# find the tfel library
if(TFEL_INSTALL_PATH)
  set(KRATOS_TFELHOME "${TFEL_INSTALL_PATH}")
else(TFEL_INSTALL_PATH)
  set(KRATOS_TFELHOME $ENV{TFELHOME})
endif(TFEL_INSTALL_PATH)

find_program(KRATOS_MFRONT       mfront
  HINTS "${KRATOS_TFELHOME}/bin")
find_program(KRATOS_TFEL_CHECK   tfel-check
  HINTS "${KRATOS_TFELHOME}/bin")
find_program(KRATOS_TFEL_CONFIG  tfel-config
  HINTS "${KRATOS_TFELHOME}/bin")

IF(NOT (KRATOS_TFEL_CONFIG AND KRATOS_MFRONT))
  MESSAGE(STATUS "tfel has not been found. Tests of the MGISApplication are disabled. "
    "Either specify the `TFEL_INSTALL_PATH` variable at the `cmake` invokation, define "
    "the `TFELHOME` variable environment or update the `PATH` and/or `LD_LIBRARY_PATH` "
    "environment variable.")
ELSE(NOT (KRATOS_TFEL_CONFIG AND KRATOS_MFRONT))
  EXECUTE_PROCESS(COMMAND ${KRATOS_TFEL_CONFIG} "--includes"
    OUTPUT_VARIABLE KRATOS_TFEL_INCLUDE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  STRING(LENGTH ${KRATOS_TFEL_INCLUDE_PATH}  KRATOS_TFEL_INCLUDE_PATH_LENGTH)
  MATH(EXPR KRATOS_TFEL_INCLUDE_PATH_LENGTH "${KRATOS_TFEL_INCLUDE_PATH_LENGTH} - 2")
  STRING(SUBSTRING ${KRATOS_TFEL_INCLUDE_PATH} 2 ${KRATOS_TFEL_INCLUDE_PATH_LENGTH}
    KRATOS_TFEL_INCLUDE_PATH)
  EXECUTE_PROCESS(COMMAND ${KRATOS_TFEL_CONFIG} "--libs"
    OUTPUT_VARIABLE KRATOS_TFEL_LIBRARY_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  STRING(LENGTH ${KRATOS_TFEL_LIBRARY_PATH}  KRATOS_TFEL_LIBRARY_PATH_LENGTH)
  MATH(EXPR KRATOS_TFEL_LIBRARY_PATH_LENGTH "${KRATOS_TFEL_LIBRARY_PATH_LENGTH} - 2")
  STRING(SUBSTRING ${KRATOS_TFEL_LIBRARY_PATH} 2 ${KRATOS_TFEL_LIBRARY_PATH_LENGTH}
    KRATOS_TFEL_LIBRARY_PATH)

  EXECUTE_PROCESS(COMMAND ${KRATOS_TFEL_CONFIG} "--cxx-standard"
    RESULT_VARIABLE KRATOS_TFEL_CXX_STANDARD_AVAILABLE
    OUTPUT_VARIABLE KRATOS_TFEL_CXX_STANDARD
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT KRATOS_TFEL_CXX_STANDARD_AVAILABLE EQUAL 0)
    set(KRATOS_TFEL_CXX_STANDARD 11)
  endif(NOT KRATOS_TFEL_CXX_STANDARD_AVAILABLE EQUAL 0)

  EXECUTE_PROCESS(COMMAND ${KRATOS_TFEL_CONFIG} "--compiler-flags"
    OUTPUT_VARIABLE KRATOS_TFEL_COMPILER_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(KRATOS_TFEL_COMPILER_FLAGS)
  
  EXECUTE_PROCESS(COMMAND ${KRATOS_TFEL_CONFIG} "--oflags"
    OUTPUT_VARIABLE KRATOS_TFEL_OFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(KRATOS_TFEL_OFLAGS)
  
  macro(find_tfel_library name)
    find_library("KRATOS_${name}"
      NAMES ${name}
      HINTS ${KRATOS_TFEL_LIBRARY_PATH})
    if(NOT KRATOS_${name})
      MESSAGE(FATAL_ERROR "${name} library not found")
    endif(NOT KRATOS_${name})
  endmacro(find_tfel_library name)

  find_tfel_library(TFELTests)
  find_tfel_library(TFELException)
  find_tfel_library(TFELUtilities)
  find_tfel_library(TFELMath)
  find_tfel_library(TFELMaterial)
  find_tfel_library(TFELPhysicalConstants)

  MESSAGE(STATUS "tfel-config           : ${KRATOS_TFEL_CONFIG}")
  MESSAGE(STATUS "mfront                : ${KRATOS_MFRONT}")
  MESSAGE(STATUS "tfel C++ standard     : ${KRATOS_TFEL_CXX_STANDARD}")
  MESSAGE(STATUS "tfel compiler flags   : ${KRATOS_TFEL_COMPILER_FLAGS}")
  MESSAGE(STATUS "tfel optimised flags  : ${KRATOS_TFEL_OFLAGS}")
  if(KRATOS_TFEL_CHECK)
    MESSAGE(STATUS "tfel-check            : ${KRATOS_TFEL_CHECK}")
  endif(KRATOS_TFEL_CHECK)  
  MESSAGE(STATUS "tfel include          : ${KRATOS_TFEL_INCLUDE_PATH}")
  MESSAGE(STATUS "tfel libs             : ${KRATOS_TFEL_LIBRARY_PATH}")
  MESSAGE(STATUS "TFELTests             : ${KRATOS_TFELTests}")
  MESSAGE(STATUS "TFELTests             : ${KRATOS_TFELTests}")
  MESSAGE(STATUS "TFELException         : ${KRATOS_TFELException}")
  MESSAGE(STATUS "TFELUtilities         : ${KRATOS_TFELUtilities}")
  MESSAGE(STATUS "TFELMath              : ${KRATOS_TFELMath}")	
  MESSAGE(STATUS "TFELMaterial          : ${KRATOS_TFELMaterial}")
  MESSAGE(STATUS "TFELPhysicalConstants : ${KRATOS_TFELPhysicalConstants}") 
  SET(HAVE_TFEL ON)

  function(kratos_add_mfront_behaviour_sources lib mat file)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront.in")
      set(mfront_file   "${CMAKE_CURRENT_BINARY_DIR}/${file}.mfront")
      configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront.in"
	"${CMAKE_CURRENT_BINARY_DIR}/${file}.mfront"
	@ONLY)
    else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront.in")
      set(mfront_file   "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront")
    endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront.in")
    set(mfront_output1 "generic/src/${file}.cxx")
    set(mfront_output2 "generic/src/${file}-generic.cxx")
    add_custom_command(
      OUTPUT  "${mfront_output1}"
      OUTPUT  "${mfront_output2}"
      COMMAND "${KRATOS_MFRONT}"
      ARGS    "--interface=generic"
      ARGS     "${mfront_file}"
      DEPENDS "${mfront_file}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generic"
      COMMENT "mfront source ${mfront_file} for interface generic")
    set(${lib}_SOURCES ${mfront_output1} ${mfront_output2}
      ${${lib}_SOURCES} PARENT_SCOPE)
  endfunction(kratos_add_mfront_behaviour_sources)

  function(kratos_mfront_behaviours_library_internal mat)
    set ( _CMD SOURCES )
    set ( _SOURCES )
    if((KRATOS_TFEL_CXX_STANDARD GREATER 17) OR (KRATOS_TFEL_CXX_STANDARD EQUAL 17))
      set(KRATOS_TFEL_MFRONT_LIBRARIES
	"${KRATOS_TFELException};${KRATOS_TFELMath};${KRATOS_TFELMaterial};${KRATOS_TFELUtilities}")
    else((KRATOS_TFEL_CXX_STANDARD GREATER 17) OR (KRATOS_TFEL_CXX_STANDARD EQUAL 17))
      set(KRATOS_TFEL_MFRONT_LIBRARIES
	"${KRATOS_TFELException};${KRATOS_TFELMath};${KRATOS_TFELMaterial};${KRATOS_TFELUtilities};${KRATOS_TFELPhysicalConstants}")
    endif((KRATOS_TFEL_CXX_STANDARD GREATER 17) OR (KRATOS_TFEL_CXX_STANDARD EQUAL 17))
    foreach ( _ARG ${ARGN})
      if ( ${_ARG} MATCHES SOURCES )
	set ( _CMD SOURCES )
      else ()
	if ( ${_CMD} MATCHES SOURCES )
          list ( APPEND _SOURCES "${_ARG}" )
	endif ()
      endif ()
    endforeach ()
    list(LENGTH _SOURCES _SOURCES_LENGTH )
    if(${_SOURCES_LENGTH} LESS 1)
      message(FATAL_ERROR "kratos_mfront_behaviours_library : no source specified")
    endif(${_SOURCES_LENGTH} LESS 1)
    include_directories("${KRATOS_TFEL_INCLUDE_PATH}")
    set(lib "${mat}Behaviours-generic")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generic")
    foreach(source ${_SOURCES})
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
	set(${lib}_SOURCES ${source} ${${lib}_SOURCES})
      else(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
	kratos_add_mfront_behaviour_sources(${lib} ${mat} ${source})
      endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    endforeach(source ${_SOURCES})
    list(LENGTH ${lib}_SOURCES nb_sources)
    if(nb_sources GREATER 0)
      message(STATUS "Adding library : ${lib} (${${lib}_SOURCES})")
      add_library(${lib} SHARED ${${lib}_SOURCES})
      target_include_directories(${lib}
	PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generic/include"
	PRIVATE "${KRATOS_TFEL_INCLUDE_PATH}")
      target_compile_options(${lib}
	PRIVATE ${KRATOS_TFEL_COMPILER_FLAGS} ${KRATOS_TFEL_OFLAGS})
      if(WIN32)
	if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	  set_target_properties(${lib}
	    PROPERTIES LINK_FLAGS "-Wl,--kill-at -Wl,--no-undefined")
	endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	install(TARGETS ${lib} DESTINATION bin)
      else(WIN32)
	install(TARGETS ${lib} DESTINATION lib${LIB_SUFFIX})
      endif(WIN32)
      target_link_libraries(${lib} ${KRATOS_TFEL_MFRONT_LIBRARIES})
    else(nb_sources GREATER 0)
      message(STATUS "No sources selected for "
	"library ${lib} for interface generic")
    endif(nb_sources GREATER 0)

    # > output library name
    set(KRATOS_MFRONT_LIB ${lib} PARENT_SCOPE)
    # <
  endfunction(kratos_mfront_behaviours_library_internal)
ENDIF(NOT (KRATOS_TFEL_CONFIG AND KRATOS_MFRONT))
  
function(kratos_mfront_behaviours_library mat)
  IF(KRATOS_TFEL_CONFIG AND KRATOS_MFRONT AND KRATOS_BUILD_TESTING)
    kratos_mfront_behaviours_library_internal(${mat} ${ARGN})
    # set_target_properties(${KRATOS_MFRONT_LIB}
    #   PROPERTIES COMPILE_OPTIONS "${KRATOS_COMPILE_OPTIONS}")
  ENDIF(KRATOS_TFEL_CONFIG AND KRATOS_MFRONT AND KRATOS_BUILD_TESTING)
endfunction()

# Install script for directory: /home/gabriela/Software/Kratos/applications/TrilinosApplication

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so"
         RPATH "~/Software/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/libs" TYPE SHARED_LIBRARY FILES "/home/gabriela/Software/Kratos/cmake_build/applications/TrilinosApplication/libKratosTrilinosCore.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so"
         OLD_RPATH "/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1:/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1:/usr/lib/x86_64-linux-gnu/libtrilinos_epetra.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscomm.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscore.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosnumerics.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosparameterlist.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosremainder.so:/usr/lib/x86_64-linux-gnu/libtrilinos_triutils.so:/usr/lib/x86_64-linux-gnu/libtrilinos_amesos.so:/usr/lib/x86_64-linux-gnu/libtrilinos_aztecoo.so:/usr/lib/x86_64-linux-gnu/libtrilinos_ifpack.so:/usr/lib/x86_64-linux-gnu/libtrilinos_loca.so:/usr/lib/x86_64-linux-gnu/libtrilinos_ml.so:/usr/lib/x86_64-linux-gnu/libtrilinos_nox.so:/usr/lib/x86_64-linux-gnu/libtrilinos_noxepetra.so:/usr/lib/x86_64-linux-gnu/libtrilinos_epetraext.so:/usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.so:/home/gabriela/Software/Kratos/cmake_build/kratos/mpi:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/gabriela/Software/Kratos/cmake_build/kratos:"
         NEW_RPATH "~/Software/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/libKratosTrilinosCore.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so"
         RPATH "~/Software/Kratos/libs")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/libs" TYPE MODULE FILES "/home/gabriela/Software/Kratos/cmake_build/applications/TrilinosApplication/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so"
         OLD_RPATH "/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1:/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1:/usr/lib/x86_64-linux-gnu/libtrilinos_epetra.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscomm.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchoscore.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosnumerics.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosparameterlist.so:/usr/lib/x86_64-linux-gnu/libtrilinos_teuchosremainder.so:/usr/lib/x86_64-linux-gnu/libtrilinos_triutils.so:/usr/lib/x86_64-linux-gnu/libtrilinos_amesos.so:/usr/lib/x86_64-linux-gnu/libtrilinos_aztecoo.so:/usr/lib/x86_64-linux-gnu/libtrilinos_ifpack.so:/usr/lib/x86_64-linux-gnu/libtrilinos_loca.so:/usr/lib/x86_64-linux-gnu/libtrilinos_ml.so:/usr/lib/x86_64-linux-gnu/libtrilinos_nox.so:/usr/lib/x86_64-linux-gnu/libtrilinos_noxepetra.so:/usr/lib/x86_64-linux-gnu/libtrilinos_epetraext.so:/usr/lib/x86_64-linux-gnu/libtrilinos_zoltan.so:/home/gabriela/Software/Kratos/cmake_build/applications/TrilinosApplication:"
         NEW_RPATH "~/Software/Kratos/libs")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/libs/KratosTrilinosApplication.cpython-36m-x86_64-linux-gnu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/KratosMultiphysics/TrilinosApplication" TYPE FILE RENAME "__init__.py" FILES "/home/gabriela/Software/Kratos/applications/TrilinosApplication/TrilinosApplication.py")
endif()


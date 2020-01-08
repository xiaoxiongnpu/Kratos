#!/bin/sh

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
# to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
# so we should add here only machine specific stuff

# an effort is made to autodetect all of the libraries needed HOWEVER:
# METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
# TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

clear

# you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                                                                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                                     \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                                   \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 "                           \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 "                                          \
-DBOOST_ROOT="~/Software/boost_1_70_0"                                              \
-DCMAKE_INSTALL_RPATH="~/Software/Kratos/libs"                                      \
-DPYTHON_EXECUTABLE="/usr/bin/python3"                                              \
-DCMAKE_BUILD_TYPE=Release                                                        \
-DMAPPING_APPLICATION=ON                                                            \
-DEXTERNAL_SOLVERS_APPLICATION=ON                                                   \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON                                               \
-DFLUID_DYNAMICS_APPLICATION=ON                                                     \
-DMESH_MOVING_APPLICATION=ON                                                        \
-DCO_SIMULATION_APPLICATION=ON                                                        \
-DMETIS_APPLICATION=ON                                                              \
-DUSE_METIS_5=ON                                                                    \
-DPARMETIS_ROOT_DIR="/home/youruser/compiled_libraries/ParMetis-3.1.1"              \
-DTRILINOS_APPLICATION=ON                                                           \
-DTRILINOS_ROOT="~/Software/trilinos"                                               \
-DTRILINOS_INCLUDE_DIR="/usr/include/trilinos"                                      \
-DTRILINOS_LIBRARY_DIR="/usr/lib/x86_64-linux-gnu"                                  \
-DTRILINOS_LIBRARY_PREFIX="trilinos_"                                               \
-DINSTALL_EMBEDDED_PYTHON=ON
-MPI_NEEDED=ON

# decomment this to have it verbose
# make VERBOSE=1 -j4
make -j8
make install

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
#to use it, copy it to another file and run it

clear

#you may want to decomment this the first time you compile
rm CMakeCache.txt

cmake ..  					\
-DCMAKE_C_COMPILER=gcc				\
-DCMAKE_CXX_COMPILER=g++			\
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -funroll-loops -ffast-math -msse3 " \
-DCMAKE_CXX_FLAGS="${CMAKE_C_FLAGS} -funroll-loops -ffast-math -msse3 " \
-DCMAKE_BUILD_TYPE=Release  			\
-DINCOMPRESSIBLE_FLUID_APPLICATION=ON  		\
-DMESHING_APPLICATION=ON 			\
-DEXTERNAL_SOLVERS_APPLICATION=ON		\
-DPFEM_APPLICATION=ON 				\
-DSTRUCTURAL_APPLICATION=ON 			\
-DCONVECTION_DIFFUSION_APPLICATION=ON 		\
-DFLUID_DYNAMICS_APPLICATION=ON 		\
-DALE_APPLICATION=ON 				\
-DFSI_APPLICATION=ON 

#decomment this to have it verbose
make VERBOSE=1 -j4
# make -j4
make install

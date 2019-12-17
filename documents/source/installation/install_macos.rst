Installation MacOS
==================

== Introduction ==

This section is devoted to explain the necessary steps for the installation of Kratos in macOS. The tutorial has been tested in all operating systems from OS X El Capitan to macOS Catalina.

We assume that the user has a basic knowledge of macOS and Unix (navigate between directories in Terminal, copy files, execute commands...).

== Install Xcode ==

The first step is to obtain the latest version of Xcode. Open the App Store and look for "xcode" in the search box. Select Xcode and start the installation. Bear in mind that all the packages require almost 6 GB of disk space and so this step can take quite a lot of time. 

If you had previously installed Xcode, it is advisable to open the App Store and go to Updates to check whether you have the latest version of Xcode.

Once the installation finishes, open Xcode in order to agree the terms of the software license agreements and complete the installation of components. After that you can quit Xcode (''cmd+q''). Before going on, it is also recommended to verify that the Command Line Tools are properly installed. To do so, open a Terminal (press ''cmd+space'', write "terminal" and hit ''enter'') and use the following command:

 xcode-select --install

If a window appears asking whether you want to install the tools now, click on "Install". You must agree the Command Line Tools License Agreement to proceed with the installation. When Xcode is properly installed, you will see the following message every time you run the "xcode-select --install" command.

 xcode-select: error: command line tools are already installed, use "Software Update" to install updates

== Install Clang with OMP support ==

Since the default Clang compiler on Xcode does not have an OpenMP implementation, here we will install a Clang version compatible with OMP. First, you need to install Homebrew (http://brew.sh), a package manager for macOS. Paste the following line at a Terminal and follow the instructions:

 /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

Once the installation is completed, type this command:

 brew install llvm

Agree if you are asked to download necessary additional packages and, after a while, you should have LLVM installed in "/usr/local/opt/llvm" and the Clang version compatible with OMP in "/usr/local/opt/llvm/bin".

== Install Python  ==

Download the latest version of Python 3 compatible with your mac from https://www.python.org/downloads/mac-osx/ . In this tutorial we used Python 3.6.6. Right-click on the downloaded file, e.g. "python-3.6.6-macosx10.9.pkg", and open it to initialise the installation wizard of Python. Once Python 3 is installed, you can erase the downloaded package to save space. Moreover, to check that Python 3 is properly configured, open a Terminal and run the following command to see the version 

 python3 --version

== Download Boost  ==

Kratos Multiphysics needs Boost libraries to support some of its functions. Download the latest compatible version of Boost from its official website http://www.boost.org/users/history/ .

We recommend you to use boost 1.67 or newer, earlier versions may cause compiling errors.

It's very important to add the correct path to the boost library in the configure.sh by setting the variable -DBOOST_ROOT. You will see an example in the configure section.

== Install CMake  ==

CMake is the tool used to compile Kratos. Download the CMake latest release for your Mac version in https://cmake.org/download/ . Open the downloaded file, e.g. "cmake-3.12.2-Darwin-x86_64.dmg", and agree the License Agreement terms. Then drag and drop the CMake application into your Applications folder.  After that you can eject "cmake-3.12.2-Darwin-x86_64" and erase it to save space.

Finally go to your Applications folder and right-click CMake.app to open it. Then you can quit CMake (''cmd+q'').

== Download Kratos  ==

Next we must download the Kratos Multiphysics source code by means of a git manager. The Command Line Tools installed with Xcode already include git, thus we just need to open a Terminal and use the following command:

 git clone https://github.com/KratosMultiphysics/Kratos Kratos

After a while, you should have the latest revision of Kratos in the folder "/Users/''youruser''/Kratos".

== Configure Kratos ==

When compiling Kratos for the first time, one needs to properly customise the configuration file. Open a Terminal, navigate to your "Kratos/cmake_build" folder and make a copy of the template file:

 cd
 cd Kratos/cmake_build/
 cp example_configure.sh.do_not_touch configure.sh

Open "configure.sh" file with any text editor and modify it like follows.  First, indicate the full path of CMake:

 /Applications/CMake.app/Contents/bin/cmake ..                                                \

Next, specify the path of the Clang compiler installed with Homebrew:

 -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang	   	  		             \

 -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++                                         \

In order to detect OpenMP_FLAG during the compilation, it is necessary to specify the path to the libraries for the Clang compiler in the CMAKE_CXX_FLAGS and CMAKE_C_FLAGS:

 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -L/usr/local/opt/llvm/lib"           \

 -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -L/usr/local/opt/llvm/lib"                          \

Now, provide the path of Boost and Python, and the location of Lapack and Blas libraries. If you followed the steps in this tutorial, just add the following lines '''setting the path to Boost properly''':

 -DBOOST_ROOT="/full/path/to/boost_1_67_0"                        \
 -DPYTHON_EXECUTABLE="/Library/Frameworks/Python.framework/Versions/3.6/bin/python3"          \
 -DLAPACK_LIBRARIES="/usr/lib/liblapack.dylib"                                                \
 -DBLAS_LIBRARIES="/usr/lib/libblas.dylib"                                                    \

Note that you can also specify which applications are going to be compiled by setting them to ON or OFF.

In the end, the configuration file could look like this:

 /Applications/CMake.app/Contents/bin/cmake ..                                                \
 -DCMAKE_C_COMPILER=/usr/local/opt/llvm/bin/clang	   	  		              \
 -DCMAKE_INSTALL_RPATH="/Users/youruser/Kratos/libs"                                          \
 -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                                     \
 -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++                                         \
 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -L/usr/local/opt/llvm/lib"           \
 -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -L/usr/local/opt/llvm/lib"                          \
 -DCMAKE_BUILD_TYPE=Release  							              \
 -DBOOST_ROOT="/full/path/to/boost_1_67_0"                                                    \
 -DPYTHON_EXECUTABLE="/Library/Frameworks/Python.framework/Versions/3.6/bin/python3"          \
 -DLAPACK_LIBRARIES="/usr/lib/liblapack.dylib"                                                \
 -DBLAS_LIBRARIES="/usr/lib/libblas.dylib"                                                    \
 -DINCOMPRESSIBLE_FLUID_APPLICATION=OFF  						      \
 -DMESHING_APPLICATION=OFF 							              \
 -DEXTERNAL_SOLVERS_APPLICATION=ON						              \
 -DPFEM_APPLICATION=OFF 							              \
 -DCONVECTION_DIFFUSION_APPLICATION=ON 						              \
 -DINSTALL_EMBEDDED_PYTHON=ON                                                                 \

'''Warning''': Cmake requires all definitions in a single line! Therefore, the line concatenation character '\' MUST NOT be followed by any whitespace in the same line as this would prevent cmake from running the lines below.

== Compile Kratos ==

In order to compile Kratos, you just need to run the following command from the same "Kratos/cmake_build/" folder:

 sh configure.sh

Please, bear in mind that depending on which applications are being compiled, the first compilation of Kratos may take a lot of time.

=== Compilation issues ===

We are aware of three warnings appearing after each linking of libraries:

 ld: warning: path '/Library/Frameworks/Python.framework/Versions/3.6/Python' following -L not a directory
 ld: warning: path '/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib' following -L not a directory
 ld: warning: path '/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib' following -L not a directory

Furthermore, the introduction of the path to the libraries of Clang in the CMAKE_CXX_FLAGS and CMAKE_C_FLAGS causes the following warning to appear:

 clang-3.9: warning: argument unused during compilation: '-L/usr/local/opt/llvm/lib'

One may suppress such warning by adding the flag '''"-Wno-unused-command-line-argument"''' to the configuration file. 

Moreover since the flag "-Wundefined-var-template" is set by default in Clang, one may see many other warnings appearing during the compilation of Kratos. To disable them, one just needs to add the flag '''"-Wno-undefined-var-template"''' to the configuration file.

Taking that into account, the CMAKE_CXX_FLAGS and CMAKE_C_FLAGS could be finally defined as follows:

 -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 -L/usr/local/opt/llvm/lib -Wno-unused-command-line-argument -Wno-undefined-var-template" \
 -DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 -L/usr/local/opt/llvm/lib -Wno-unused-command-line-argument -Wno-undefined-var-template"                \

== Set up your shell environment ==

After compiling Kratos, you need to tell macOS where to find the libraries by adding some paths to your shell environment. With macOS Catalina, execute the following commands to create or open (if it already exists) the ".zprofile":

 cd
 nano .zprofile

'''Note''': For previous macOS versions the configuration file is called ".bash_profile" instead of ".zprofile".

Now add the following two lines containing the path to your Kratos folder and the path to the compiled libraries. If you followed the steps in this tutorial they are:

 export PYTHONPATH="/Users/youruser/Kratos:$PYTHONPATH"
 export DYLD_LIBRARY_PATH="/Users/youruser/Kratos/libs:$DYLD_LIBRARY_PATH"

Save the file by pressing ''ctrl+o'' and then hitting ''enter'', exit nano with ''ctrl+x'', and quit the Terminal with ''cmd+q''. After this, every time you open a new Terminal window, the paths will be set automatically.

== Test Kratos ==

The easiest way to test if Kratos is properly installed is to reproduce the Kratos Multiphysics message. Open a Terminal and execute Python 3 with:

 python3

Once you are in Python 3, execute the following line:

 import KratosMultiphysics

If everything is properly installed you should see this message:

  |  /           |             
  ' /   __| _` | __|  _ \   __|
  . \  |   (   | |   (   |\__ \ 
 _|\_\_|  \__,_|\__|\___/ ____/
            Multi-Physics 5.0.17106

To quit Python 3 just run

 exit()

And that's all. You can now run your own python scripts to test that everything works.

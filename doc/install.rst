===============================
Installation for Ubuntu/Debian
===============================

Prerequisites:
 
	+ C++17 compatible compiler

	+ CMake >= 3.18
	
	+ pytorch

	+ boost graph library

	Optional - Documentation:

	+ Doxygen >= 1.9.0

	+ Sphinx >= 3.2.1 (You can install it with `pip`. Since it is python3 you may need to install `pip3`, instead of `pip`)

	+ Breathe >= 4.20.0 (You can install it with `pip` as well.)

	Optional - Python bindings:

	+ pybind11 

	+ nvidia-cuda-toolkit (May be required when developing in c++.)

	
1. Obtaining TORCHAMI:
 
	Clone Git repository:

	::

	$ git clone https://github.com/mdburke11/torchami.git

2. Obtaining Dependencies

	Required

	pytorch:

	::

	$ pip install torch 

	boost graph:

	::

	$ sudo apt-get install libboost-graph-dev

	Optional - but recommended

	::

	$ sudo apt-get install pybind11 
	$ sudo apt-get install nvidia-cuda-toolkit 

3. Building:

	An example cmake command is in the `compile.sh` script.  Open the file in any editor, and change relevant file paths, in particular -DCMAKE_INSTALL_PREFIX=/path/to/install should be modified to your chosen install path.
	Use a standard `CMake` procedure: 

			::

			$ mkdir build install
			$ cd build
			$ sh ../compile.sh 
			$ make install 



------------------------
Using with your projects
------------------------

Compilation will produce two libraries libtamigraph and libtorchami.  Calculations only require libtorchami but labelling Feynman diagrams is functionality within libtamigraph. 
Including the headers and linking to the libraries is required.  This can be accomplished using a standard `CMake`-based approach:

		::

		 
		  $ export torchami_DIR=/install/dir/of/torchami
		  $ cd /your/project/
		  $ cat >CMakeLists.txt
		  cmake_minimum_required(VERSION 3.18)
		  project(MyProject C CXX)
		  find_package(torchami REQUIRED)
		  add_executable(my_program main.cpp)
		  target_link_libraries(my_program ${torchami})
		  ...cont
		  $ cmake .
		  $ make


-------------------------
Python
-------------------------

A python library is generated using pybind11.  The installation path for this library will likely throw an error.  
After compilation with the `-DMAKE_PYTAMI=ON` flag the pytami library will be located in `/build/src/pytami_src/pytami.cpython-310-x86_64-linux-gnu.so` and can be manually installed or imported as necessary.

-------------------------
Troubleshooting
-------------------------

| Older GPUs
| (note: cuda 10.2 or greater is required, so ubuntu 22.04 or later recommended.  Otherwise install from nvidia instructions: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#ubuntu )\
| Caffe2 requires cuDNN https://developer.nvidia.com/rdp/cudnn-download

| Download keys from above then:
| sudo apt-get install libcudnn8=8.9.2.26-1+cuda12.1            (or appropriate versions)
| sudo apt-get install libcudnn8-dev=8.9.2.26-1+cuda12.1

| May need: -DCMAKE\_CUDA\_COMPILER:PATH=/usr/local/cuda/bin/nvcc

Then just make with the compile script as above

------------

Note: When compiling executables in c++, if the python torch libraries are installed via python's pip command they may be installed in a non-standard location such as: /home/username/.local/lib/python3.10/site-packages/torch/lib/ . It this is the case then you can just add this directory to LD\_LIBRARY\_PATH.  Same holds for pybind11. 

------------

IMPORTANT Note: If using pytorch from the pip wheel, linking a c++ code can be a bit confusing.  
pytorch is compiled with a c++ flag \_GLIBCXX\_USE\_CXX11\_ABI=0 . 
Any third party library that is not compiled with this flag will be ignored by the linker.  
Solution for c++ is to compile all your DIRECTLY linked libraries with this flag. 


===============================
Installation for MACOS
===============================

Since MACOS does not support cuda, we do not recommend using torchami.  
Nevertheless, the pytorch dependency can be met on MACOS via the cpu-only build of pytorch.

| pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

We do not recommend using `torchami` on non-gpu systems. For CPU only cases the libami code is preferred: https://github.com/jpfleblanc/libami .
	
.. _`Github wiki`: https://github.com/mdburke11/torchami

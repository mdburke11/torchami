============
Installation
============

For detailed instructions, please see `Github Wiki`_

Prerequisites:
 
	+ C++17 compatible compiler

	+ CMake >= 3.18

	Optional:

	+ Doxygen >= 1.9.0

	+ Sphinx >= 3.2.1 (You can install it with `pip`. Since it is python3 you may need to install `pip3`, instead of `pip`)

	+ Breathe >= 4.20.0 (You can install it with `pip` as well.)

	+ pybind11 

	+ nvidia-cuda-toolkit 

	+ pytorch

Useful commands:\
sudo apt-get install pybind11 \
sudo apt-get install nvidia-cuda-toolkit 

Python 3:\
pip install torch  \
(note: if python isn't found try: sudo apt-get install python-is-python3 )


1. Obtaining TORCHAMI:
 
	Clone Git repository:

	::

	$ git clone https://github.com/mdburke11/torchami.git

2. Building:

	An example cmake command is in the `compile.sh` script.  Open the file in any editor, and change relevant file paths.
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




	
.. _`Github wiki`: https://github.com/mdburke11/torchami

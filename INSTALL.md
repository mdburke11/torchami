
Installation
------------

Things you'll need:\
sudo apt-get install pybind11 \
sudo apt-get install nvidia-cuda-toolkit 

Python 3:\
pip install torch  \
(note: if python isn't found try: sudo apt-get install python-is-python3 )

Compilation:\
cd torchami\
mkdir build install\
cd build \
sh ../compile.sh \
make test \
make install 

------------

Older GPUs\
(note: cuda 10.2 or greater is required, so ubuntu 22.04 or later recommended.  Otherwise install from nvidia instructions: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#ubuntu )\
Caffe2 requires cuDNN https://developer.nvidia.com/rdp/cudnn-download

Download keys from above then:\
sudo apt-get install libcudnn8=8.9.2.26-1+cuda12.1            (or appropriate versions)\
sudo apt-get install libcudnn8-dev=8.9.2.26-1+cuda12.1

May need: -DCMAKE\_CUDA\_COMPILER:PATH=/usr/local/cuda/bin/nvcc\

Then just make with the compile script as above\

------------

Note: When compiling executables in c++, the python torch libraries may be installed in a non-standard location such as: /home/username/.local/lib/python3.10/site-packages/torch/lib/ . It this is the case then you can just add this directory to LD\_LIBRARY\_PATH.

------------

Note: If using pytorch from the pip wheel, linking a c++ code can be a bit confusing.  pytorch is compiled with a c++ flag \_GLIBCXX\_USE\_CXX11\_ABI=0 . Any third party library that is not compiled with this flag will be ignored by the linker.  Solution for c++ is to compile pytorch for yourself, or compile all your DIRECTLY linked libraries with this flag. \
	
	
.. _`Github wiki`: https://github.com/mdburke11/torchami


Installation
------------

Things you'll need:
sudo apt-get install pybind11
sudo apt-get install nvidia-cuda-toolkit

Python 3
pip install torch  
(note: if python isn't found try: sudo apt-get install python-is-python3 )

Compilation
cd torchami
mkdir build install
cd build
sh ../compile.sh 
make test
make install 

------------

Older GPUs
(note: cuda 10.2 or greater is required, so ubuntu 22.04 or later recommended.  Otherwise install from nvidia instructions: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#ubuntu )
Caffe2 requires cuDNN https://developer.nvidia.com/rdp/cudnn-download

Download keys from above then:
sudo apt-get install libcudnn8=8.9.2.26-1+cuda12.1            (or appropriate versions)
sudo apt-get install libcudnn8-dev=8.9.2.26-1+cuda12.1

May need: -DCMAKE_CUDA_COMPILER:PATH=/usr/local/cuda/bin/nvcc

Then just make with the compile script as above
	
	
.. _`Github wiki`: https://github.com/mdburke11/torchami


Installation
------------

Things you'll need:
sudo apt-get install pybind11
sudo apt-get install nvidia-cuda-toolkit
(note: cuda 10.2 or greater is required, so ubuntu 22.04 or later recommended.  Otherwise install from nvidia instructions: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#ubuntu )
Caffe2 requires cuDNN https://developer.nvidia.com/rdp/cudnn-download

Download keys from above then:
sudo apt-get install libcudnn8=8.9.2.26-1+cuda12.1            (or appropriate versions)
sudo apt-get install libcudnn8-dev=8.9.2.26-1+cuda12.1


May need: -DCMAKE_CUDA_COMPILER:PATH=/usr/local/cuda/bin/nvcc
          

Python 3

pip install torch 


Working build from Ubuntu 22.04 LTS is just

pip install torch
sudo apt-get install python-is-python3
sudo apt-get install nvidia-cuda-toolkit

Then just make with the compile script
	
	
.. _`Github wiki`: https://github.com/mdburke11/torchami

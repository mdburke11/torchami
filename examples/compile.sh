CC=$(which mpicc) CXX=$(which mpicxx) cmake -DPYTHON_EXECUTABLE=`which python` -DCMAKE_PREFIX_PATH=`python -c "import torch; print(torch.utils.cmake_prefix_path)"` ..

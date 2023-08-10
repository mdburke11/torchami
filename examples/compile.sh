CC=$(which mpicc) CXX=$(which mpicxx) cmake -DPYTHON_EXECUTABLE=`which python3` -DCMAKE_PREFIX_PATH=`python3 -c "import torch; print(torch.utils.cmake_prefix_path)"` ..

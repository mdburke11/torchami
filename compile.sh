cmake -DCMAKE_INSTALL_PREFIX=/path/to/torchami/install -DPYTHON_EXECUTABLE=`which python` -DCMAKE_PREFIX_PATH=`python -c "import torch; print(torch.utils.cmake_prefix_path)"` ..

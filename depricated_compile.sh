cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=/path/to/install -DCMAKE_BUILD_TYPE=Release -DBUILD_DOC=OFF  -DMAKE_PYTAMI=ON -DPYTHON_LIBRARY_DIR=`python -c "import site; print(site.getsitepackages()[0])"` -DPYTHON_EXECUTABLE=`which python` -DCMAKE_PREFIX_PATH=`python -c "import torch; print(torch.utils.cmake_prefix_path)"` ..


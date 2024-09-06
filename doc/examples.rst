=============
Examples
=============

Examples are contained in the ``/torchami/examples`` folder.  They can be compiled and run using the same procedure as for the main code.  The relative path to /torchami/install is set by default.  Examples are explained in detail in the code paper and contain additional comments and timing information that may be helpful.  Compilation of the examples follows:  


		::
		
		 $ cd examples
		 $ mkdir build && cd build
		 $ sh ../compile.sh
		 $ make
		 $ ./examples

There are also examples of using the python bindings, contained in ``/torchami/pytami_examples``. After compliling the library with ``MAKE_PYTAMI=ON`` these example scripts can be ran using the python environment where `pytami` was compiled.


		::
		
		 $ python example_main.py
		 $ python graph_example.py
		 $ python momentum_main.py

Further information and updates will be posted on the `Github Wiki`_. 
	
.. _`Github wiki`: https://github.com/mdburke11/torchami

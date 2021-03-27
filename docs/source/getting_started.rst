Getting Started
===============

PyPI build
^^^^^^^^^^
If you used PyPI to build the project, Download the input template to your working directory

.. code-block:: shell

  wget https://raw.githubusercontent.com/sayanadhikari/PPDyn/main/input.ini

Now, either create a python script in your working directory or use your python console

.. code-block:: python

  from PPDyn import ppdyn
  from PPDyn.ppdplot import animate
  import time

  start = time.time()
  ppdyn(input)
  end = time.time()
  print("Elapsed (after compilation) = %s"%(end - start)+" seconds")
  animate()

GNU Make build
^^^^^^^^^^^^^^
If you used GNU Make to build the project, upon successful compilation, run the code using following command

.. code-block:: shell

  ppdyn --i input.ini

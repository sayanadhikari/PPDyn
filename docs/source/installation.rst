Installation
============
Prerequisites
-------------
- .. _GNU Make: https://www.gnu.org/software/make/
- Python3 or higher: https://www.python.org/download/releases/3.0/
- git: https://git-scm.com/

Procedure
---------
Install the *PPDyn* using either PyPI or directly from GitHub

Using PyPI
----------
.. code-block:: shell

  pip install PPDyn

Getting Started
---------------
Download the input template to your working directory

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

Using GNU Make
--------------
First make a clone of the master branch using the following command

.. code-block:: shell

  git clone https://github.com/sayanadhikari/PPDyn.git

Then enter inside the $project directory

.. code-block:: shell

  cd PPDyn
Now complile and built the $project code

.. code-block:: shell

  make all

Getting Started
---------------
Upon successful compilation, run the code using following command

.. code-block:: shell

  ppdyn --i input.ini

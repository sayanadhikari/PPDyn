Installation
============

*PPDyn* can be installed in two ways.
#. Using *pip* from PyPI
#. Using *git clone* from GitHub

Prerequisites
-------------
- `Python3 or higher <https://www.python.org/download/releases/3.0/>`_
- `GNU Make <https://www.gnu.org/software/make/>`_ (If you are installing from GitHub)
- `git <https://git-scm.com/>`_ (If you are installing from GitHub)

Procedure
---------

.. Using *pip* from PyPI::
^^^^^^^^^^
.. code-block:: shell

  pip install PPDyn

.. Using *git clone* from GitHub::
^^^^^^^^^^^^^^
First make a clone of the master branch using the following command

.. code-block:: shell

  git clone https://github.com/sayanadhikari/PPDyn.git

Then enter inside the *PPDyn* directory

.. code-block:: shell

  cd PPDyn

Now complile and built the *PPDyn* code

.. code-block:: shell

  make all

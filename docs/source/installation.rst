Installation
============
Prerequisites
-------------
- `Python3 or higher <https://www.python.org/download/releases/3.0/>`_
- `GNU Make <https://www.gnu.org/software/make/>`_ (If you are installing from GitHub)
- `git <https://git-scm.com/>`_ (If you are installing from GitHub)

Procedure
---------
*PPDyn* can be installed in two ways.

#. `Using *pip* from PyPI <#id1>`_
#. `Using *git clone* from GitHub <#id2>`_

Using *pip* from PyPI
^^^^^^^^^^^^^^^^^^^^^
.. code-block:: bash

  pip install PPDyn

Using *git clone* from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First make a clone of the master branch using the following command

.. code-block:: bash

  git clone https://github.com/sayanadhikari/PPDyn.git

Then enter inside the *PPDyn* directory

.. code-block:: bash

  cd PPDyn

Now complile and built the *PPDyn* code

.. code-block:: bash

  make all

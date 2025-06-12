.. _installation:

Installation
============

.. contents:: :local:

.. warning::

   The following is a set of instructions to get Noda up and running
   *in most cases*. Issues may arise depending on how your Python installation
   is configured. You may also want to follow a different installation route.
   As with any Python package, it is advisable to have a general knowledge of
   virtual environments, packages and dependencies before attempting to install
   Noda.

Installing Noda
---------------

Download and unzip the package archive. Then open a terminal and go to
the unzipped folder (the package root directory aka installation directory):

.. tab:: Unix (Linux / Mac)

   .. code-block::

      $ cd /path/to/NODA

.. tab:: Windows

   .. code-block::

      > cd \path\to\NODA

It is recommended to install Noda in a virtual environment. These instructions
assume you use a virtual environment named noda-env. It can be installed using
pip (recommended) or conda:

* with pip:

    .. tab:: Unix (Linux / Mac)

       .. code-block::

          $ python -m venv noda-env          # create the virtual environment
          $ source noda-env/bin/activate     # activate the virtual environment

    .. tab:: Windows

       .. code-block::

          > python -m venv noda-env          # create the virtual environment
          > .\noda-env\Scripts\activate      # activate the virtual environment

* with conda:

    .. tab:: Unix (Linux / Mac)

       .. code-block::

          $ conda create -n noda-env         # create the virtual environment
          $ source activate noda-env         # activate the virtual environment

    .. tab:: Windows

       .. code-block::

          > conda create -n noda-env         # create the virtual environment
          > conda activate noda-env          # activate the virtual environment

In both cases, it is recommended to install Noda with pip:

.. code-block::
   
   (noda-env) $ pip install .

.. _running_tests:

Running tests
-------------

After Noda is installed, it is recommended to run the test suite, which relies
on pytest:

.. code-block::

      (noda-env) $ pip install pytest
      (noda-env) $ pytest

Warnings may be issued and can be ignored, but all tests should pass.

Building the documentation
--------------------------

It is also recommended to build the documentation locally. This is done with
Sphinx. To install Sphinx:

.. code-block::
   
   (noda-env) $ pip install sphinx
   (noda-env) $ pip install sphinx_rtd_theme
   (noda-env) $ pip install sphinx-inline-tabs

Once you have Sphinx set up, go to the ``doc`` folder of the installation
directory and build the documentation:

.. code-block::

   (noda-env) $ cd doc
   (noda-env) $ make html
      
The entry point to the documentation, ``index.html``, is located in
``doc/_build/html``. It is useful to make a link or a shortcut to this file
and place it somewhere easier to access:

.. tab:: Unix (Linux / Mac)

   .. code-block::

      (noda-env) $ ln -s /path/to/NODA/doc/_build/html/index.html ~/Desktop/noda_documentation

.. tab:: Windows

   .. code-block::

      (noda-env) > $WshShell = New-Object -ComObject WScript.Shell
      (noda-env) > $Shortcut = $WshShell.CreateShortcut("C:\Users\user\Desktop\noda_documentation.lnk")
      (noda-env) > $Shortcut.TargetPath = "C:\path\to\NODA\doc\_build\html\index.html"
      (noda-env) > $Shortcut.Save()

.. _environment_variable:

Setting the ``NODA_HOME`` environment variable
----------------------------------------------

This is used to tell the package where user data files (such as thermodynamic
databases) will be located.

.. tab:: Unix (Linux / Mac)

   Add an export command in your Bash shell startup script (``~/.bashrc``):
   
   .. code-block::
   
      export NODA_HOME=/the/path/you/want

.. tab:: Windows

   Type "environment variables" in the search tool, then in the Environment
   variables window, choose "New", then in Variable name, type NODA_HOME, and in
   Variable value, choose the path you want.

Uninstalling Noda
-----------------

.. code-block::

   (noda-env) $ pip uninstall noda

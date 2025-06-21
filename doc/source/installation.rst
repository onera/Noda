.. _installation:

Installation
============

.. contents:: :local:

.. warning::

   The following is a set of instructions to get Noda up and running
   *in most cases*. It assumes you already have a Python installation with pip
   accessible from a command-line interface. Issues may arise depending on
   how your internet connection, Git and Python installations are configured.
   You may also want to follow a different installation route. As with any
   Python package, it is advisable to have a general knowledge of virtual
   environments, packages and dependencies before attempting to install Noda.

Getting, installing and updating Noda
-------------------------------------

It is recommended to install Noda in a virtual environment. These instructions
assume you use a virtual environment named noda-env. It can be installed using
the standard library venv module or with other tools such as conda. Open a
terminal, cd to a folder where you want to keep the installation directory,
create and activate the virtual environment:

* with venv:

    .. tab:: Unix (Linux / Mac)

       .. code-block::

          $ python -m venv noda-env
          $ source noda-env/bin/activate

    .. tab:: Windows

       .. code-block::

          > python -m venv noda-env
          > .\noda-env\Scripts\activate

* with conda:

    .. tab:: Unix (Linux / Mac)

       .. code-block::

          $ conda create -n noda-env
          $ source activate noda-env

    .. tab:: Windows

       .. code-block::

          > conda create -n noda-env
          > conda activate noda-env

Then get and install Noda via either of the following routes:

Installing from the url
^^^^^^^^^^^^^^^^^^^^^^^

Install directly from the repository url:

.. code-block::
   
   (noda-env) $ pip install git+https://github.com/onera/Noda.git@main

To update, run the same command.

Cloning the repository
^^^^^^^^^^^^^^^^^^^^^^

First clone the repository, then install:

.. code-block::
   
   (noda-env) $ git clone https://github.com/onera/Noda.git
   (noda-env) $ cd Noda
   (noda-env) $ pip install .

To update, pull the repository and run the install command.

Downloading the package archive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use a web browser to find the online repository: https://github.com/onera/noda.
Click on the green "Code" button, then "Download ZIP". Copy the archive where
you put the virtual environment, then unzip and install:

.. code-block::

   (noda-env) $ tar -xf Noda-main.zip
   (noda-env) $ cd Noda-main
   (noda-env) $ pip install .

To update, download the archive again and repeat.

.. hint::

   If the ``pip install git+https`` command fails or if the ``git clone``
   command fails with the following error:
   ``fatal: unable to access 'https://github.com/onera/Noda.git/': Failed to connect to github.com port 443 after 21064 ms: Could not connect to server``,
   this may be because you are connecting to the web via a proxy server and Git
   is not aware of it. To configure Git to use a proxy, run:

   .. code-block::

      $ git config --global http.proxy yourhost:yourport
   
   Your proxy host and port can be found in a number of ways, see for instance
   https://superuser.com/q/346372. Alternatively, you can try disconnecting
   your VPN if you were using one.

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
Sphinx. To install Sphinx and the necessary extensions:

.. code-block::
   
   (noda-env) $ pip install sphinx sphinx_rtd_theme sphinx-inline-tabs

Then go to the ``doc`` folder of the installation directory and build the
documentation:

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

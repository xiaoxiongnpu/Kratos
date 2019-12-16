Getting Kratos Binaries for Linux
=================================

This page covers the process of downloading a compiled version of Kratos for Linux and how to use it to run an example. 

Minimum requirements
--------------------
Ubuntu 14.04 or greater. Other flavors of Linux are supported but not tested.

Download and Installation
-------------------------
You can find the Kratos binaries for Linux in the `Release Section <https://github.com/KratosMultiphysics/Kratos/releases/tag/7.0>`_:

1) Download the tgz file for linux

2) Extract the contents of the tgz in a location of your choice

3) Recommended: Execute the :code:`install.sh` script

This finalizes the installation process.

Usage
-----
- Download any example of your interest from the `Examples Repository <https://github.com/KratosMultiphysics/Examples>`_ and execute:


.. code-block:: bash

   >kratos MainKratos.py


Advanced
........
If you prefer to skip the installation part, or have multiple versions of kratos installed simultaneously:

1) Copy this :code:`runner.sh` script and edit the contents of the :code:`KRATOS_ROOT` variable:

.. code-block:: bash

    KRATOS_ROOT=/path/to/kratos   # Please change this variable
    export PYTHONHOME=${KRATOS_ROOT}
    export PYTHONPATH=${KRATOS_ROOT}:${KRATOS_ROOT}/python34.zip
    export LD_LIBRARY_PATH=${KRATOS_ROOT}/libs:${KRATOS_ROOT}/OpenMPI/lib:/home/roigcarlo/KratosInstall/libs
    ${KRATOS_ROOT}/runkratos $*


Do not forget to give the script execution rights:

.. code-block:: bash

   >chmod +x runner.sh


2) Make an alias for the script

.. code-block:: bash

   >echo "alias kratos=runner.sh" >> $HOME/.bashrc


3) Open the terminal and execute the problem script using kratos:

.. code-block:: bash

   >kratos MainKratos.py



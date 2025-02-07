Installation
============




Currently, the installation is done only using the GitHub repository. Updated instructions will be added once the toolbox is added to MDAKits.

Vitual environment
-----------------


We strongly recommend the creation of a virtual environment with conda


.. code-block:: console

    conda create --name twod
    conda activate twod

.. tip::

    ``conda create`` command creates a conda environment named *twod*
    and ``conda activate`` command starts environment in your terminal.
    Whenever you close your terminal, you need to use activate command again to restart environment.



Installation from github
------------------------

Download the repository and go into the directory twdanalysis:

.. code-block:: console

    git clone https://github.com/monjegroup/twod-analysis-kit.git
    cd twod-analysis-kit

Inside the directory, install the development and documentation dependencies:

.. code-block:: console

    conda env update --name twod --file devtools/conda-envs/test_env.yaml
    conda env update --name twod --file docs/requirements.yaml

Finally, build the current version via pip command

.. code-block:: console

   pip install -e .

To check for the latest version of the repositoy, navigate to the twoanalysis directory and run:

.. code-block:: console

   git pull

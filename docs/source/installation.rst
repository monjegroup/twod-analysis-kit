Installation
============




Currently we only support the installation can be done using the github repository. This
will change once this repository is added to MDAKits.

Vitual environment
-----------------

.. code-block:: console

    conda create --name twod
    conda activate twod

.. tip::

    ``conda create`` command creates a conda environment named *twod*
    and ``conda activate`` command starts environment in your terminal.
    Whenever you close your terminal, you need to use activate command again to restart environment.

Before start, we strongly recommend the creation of a virtual environment with conda

Installation from github
------------------------

The following two lines download our repository and go into the directory twdanalysis

.. code-block:: console

    git clone https://github.com/monjegroup/twod-analysis-kit.git
    cd twod-analysis-kit

Now, once in the directory, install the development and documentation dependencies:

.. code-block:: console

    conda env update --name twod-analysis-kit --file devtools/conda-envs/test_env.yaml
    conda env update --name twod-analysis-kit --file docs/requirements.yaml

Finally, build dependencies and twodanalysis current version via pip command

.. code-block:: console

   pip install -e .

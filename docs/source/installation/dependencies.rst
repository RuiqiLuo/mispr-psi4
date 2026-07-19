===============================
Prerequisites
===============================

Conda environment
------------------------------
MISPR depends on a number of third party Python packages, and usually on
specific versions of those packages. In order not to interfere with third
party packages needed by other software on your machine or cluster, we
strongly recommend isolating MISPR in a conda environment. In the
following, we describe how to create a conda environment using
the `conda <https://docs.conda.io/projects/conda/en/latest/>`_ tool.

.. important::

   MISPR requires Python version 3.10 or higher. We have extensively tested MISPR with Python 3.10.

Creating the environment: the exact recipe
============================================
Before creating a conda environment, ensure that Anaconda or Miniconda is installed on your system.
Most HPC clusters provide Anaconda as a loadable module. If you need to install it yourself, you can
install Miniconda by following the `official installation guide <https://docs.conda.io/projects/miniconda/en/latest/>`_.

Some MISPR dependencies must come from **conda** (conda-forge) and others
from **pip**, and the **order matters** -- see the warning below for why.
Run these commands exactly in this order::

    # 1. create and activate the environment
    conda create -n mispr_env python=3.10
    conda activate mispr_env

    # 2. FIRST: everything that ships compiled binaries, from conda-forge
    #    (openbabel is always required; psi4 + dftd3-python only if you plan
    #    to use the Psi4 backend -- harmless to include either way)
    conda install -c conda-forge openbabel=3.1.1 psi4 dftd3-python

    # 3. THEN: MISPR itself via pip -- pip will see numpy/scipy/pandas already
    #    provided by conda-forge and reuse them instead of replacing them.
    #    NOTE: the PyPI "mispr" package does NOT include the ORCA/Psi4
    #    backends; install from this repository instead:
    pip install git+https://github.com/RuiqiLuo/mispr-psi4.git
    #    (or, for development mode: git clone the repository first, then
    #     pip install -e /path/to/your/clone)

    # 4. only if you will use the Psi4 ESP workflow: the resp package,
    #    which MUST come from GitHub, not from PyPI (see warning below)
    pip install git+https://github.com/cdsgroup/resp.git

After activation, your prompt should have ``(mispr_env)`` in front of it.
To deactivate the environment, simply run ``conda deactivate``.

.. warning::
   **Why the order matters -- the numpy/scipy/pandas conflict.** pip wheels
   and conda-forge builds of ``numpy``, ``scipy``, and ``pandas`` are not
   always binary-compatible with each other. If pip installs its own numpy
   first (as a dependency of mispr) and conda-forge's psi4/openbabel are
   added on top later, imports start failing with errors like::

       numpy.dtype size changed, may indicate binary incompatibility
       ValueError: All ufuncs must have type numpy.ufunc

   Installing the conda-forge packages *first* (step 2) and pip packages
   *second* (step 3) avoids this entirely. Two rules keep the environment
   healthy afterwards:

   * never run ``pip install numpy`` / ``scipy`` / ``pandas`` (or
     ``pip install --upgrade`` anything that pulls them in) in this
     environment;
   * if the environment does break with the errors above, repair it with
     ``conda install -c conda-forge numpy scipy pandas --force-reinstall``
     -- or, often faster in practice, delete the environment and redo the
     recipe from step 1 in the correct order.

.. note::
   You may need to install ``pip`` and ``setuptools`` in your conda
   environment in case the system or user version of these tools is old::

    conda install pip setuptools

Required Software
---------------------------------

.. _chem-software-deps:

Computational chemistry software dependencies
=============================================
At the backend, MISPR uses:

.. list-table:: 
   :widths: 20 25 40 15
   :header-rows: 1

   * - Software
     - License Type
     - Purpose
     - Installation
   * - `Gaussian <https://gaussian.com>`_
     - Commercial
     - Perform DFT calculations
     - License required
   * - `Psi4 <https://psicode.org>`_
     - Open Source
     - Perform DFT calculations (alternative to Gaussian)
     - ``conda install -c conda-forge psi4``
   * - `ORCA <https://www.faccts.de/orca/>`_
     - Free for academic use (registration required)
     - Perform DFT calculations (alternative to Gaussian)
     - Download from the `ORCA forum <https://orcaforum.kofo.mpg.de>`_
   * - `AmberTools24 <https://ambermd.org/AmberTools.php>`_
     - Open Source
     - Generate GAFF parameters
     - Direct download
   * - `LAMMPS <https://www.lammps.org>`_
     - Open Source
     - Run MD simulations
     - Direct download
   * - `Packmol <https://m3g.github.io/packmol/download.shtml>`_
     - Open Source
     - Build initial MD configurations
     - Follow `user guide <https://m3g.github.io/packmol/userguide.shtml>`_
   * - `Schrödinger <https://www.schrodinger.com/>`_
     - Commercial (Academic license available)
     - Generate OPLS2005 parameters (Optional)
     - License required

Ensure that you have access to the executables of these software
before using MISPR. If Gaussian, AmberTools, Schrödinger and LAMMPS are already installed on HPC
machines, the user typically needs to load their corresponding modules
before their use.

Psi4 backend (alternative to Gaussian)
=======================================
As of version 0.0.5, MISPR also supports running DFT calculations through
`Psi4 <https://psicode.org>`_ instead of Gaussian. Unlike Gaussian, Psi4 is
free/open source and is driven entirely through its Python API (no separate
executable, license, or input/output file handling required) -- so instead of
loading a module, you install it directly into your conda environment::

    conda install -c conda-forge psi4 dftd3-python

MISPR's ESP workflow additionally needs the ``resp`` package (RESP charge
fitting), which is not published on PyPI/conda-forge under that name and must
be installed from its GitHub source::

    pip install git+https://github.com/cdsgroup/resp.git

.. important::
   Do **not** ``pip install resp`` on its own -- PyPI has an unrelated,
   unofficial package also named ``resp``; installing it silently gives you
   the wrong package (no ``resp.resp()`` function), and code depending on it
   will fail with ``AttributeError: module 'resp' has no attribute 'resp'``.

.. warning::
   Install Psi4 **before** installing MISPR via pip -- i.e. follow the exact
   order given in
   :ref:`installation/dependencies:Creating the environment: the exact recipe`
   above. Installing Psi4 into an environment where pip already installed
   ``numpy``/``scipy``/``pandas`` triggers the binary-incompatibility errors
   described there.

ORCA backend (alternative to Gaussian)
=======================================
As of version 0.0.5, MISPR also supports running DFT calculations through
`ORCA <https://www.faccts.de/orca/>`_ instead of Gaussian. Like Gaussian (and
unlike Psi4), ORCA runs as an external program: MISPR writes a text input file,
invokes the ``orca`` binary, and parses the text output. ORCA is free for
academic use but requires registration: create an account on the
`ORCA forum <https://orcaforum.kofo.mpg.de>`_, download the Linux archive
matching your cluster (e.g. ``orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz``),
transfer it to your machine, and extract it::

    tar -xf orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz

.. warning::
   The extracted installation is large (~17 GB for ORCA 6.1). On HPC clusters
   with small home-directory quotas, extract it into a scratch/project
   filesystem instead of ``$HOME``.

MISPR locates the ORCA executable through (in priority order) the ``orca_cmd``
argument accepted by the ORCA workflows, the ``ORCA_CMD`` environment variable,
or a bare ``orca`` on your ``PATH``. The simplest setup is to export the full
path once::

    export ORCA_CMD=/path/to/orca_6_1_1_linux_x86-64_shared_openmpi418/orca

.. note::
   Running ORCA in parallel (``num_cores`` > 1) requires ``ORCA_CMD`` (or
   ``orca_cmd``) to be the binary's full absolute path -- an ORCA/OpenMPI
   requirement, not a MISPR one -- and the matching OpenMPI version available
   in your environment.

.. _py-package-deps:

Python package dependencies
=================================
* `pymatgen <https://pymatgen.org>`_: MISPR uses pymatgen for handling
  different molecule representations and i/o operations specific to
  Gaussian and LAMMPS. We have made changes to the pymatgen library to
  make it compatible with our needs in MISPR. These changes have not
  been merged yet with the main pymatgen library. Therefore, in order
  to use MISPR, you need to install the MolMD version of pymatgen by
  running the following commands in your ``codes`` directory::

    pip3 install pymatgen@git+https://github.com/molmd/pymatgen@molmd_fix_3-9#egg=pymatgen

* `FireWorks <https://materialsproject.github.io/fireworks/>`_: MISPR
  uses FireWorks to design, manage, and execute workflows.

  Further details can be found in the `FireWorks documentation  <https://materialsproject.github.io/fireworks/installation.html>`_.

  .. note::
   While FireWorks is used in MISPR for managing the DFT and MD
   workflows due to its many advantages, it takes some time to learn
   and get used to it.

* `custodian <https://materialsproject.github.io/custodian/>`_: MISPR uses
  custodian for handling errors that occur during the simulations and
  correcting them according to predefined rules. We have contributed a Gaussian
  plug-in to the custodian library, and these changes have been merged with 
  the main custodian library.

* `OpenBabel <https://openbabel.org>`_ to handle molecule operations 
  via pymatgen as an interface. You can install OpenBabel using conda::

    conda install -c conda-forge openbabel=3.1.1

* `MDPropTools <https://github.com/molmd/mdproptools>`_: MISPR uses mdproptools, which is a standalone 
  Python package we developed for analyzing molecular dynamics trajectories and 
  output files. 

.. note::
   FireWorks, custodian, and MDPropTools will be automatically installed as dependencies when you 
   install MISPR. You don't need to install them separately.

MongoDB
-------------------------
MISPR uses `MongoDB <https://docs.mongodb.com/manual/>`__ as the backend database.
MongoDB is a NoSQL database that is designed to store and retrieve
data in a highly efficient and scalable manner. It stores data in the
form of documents represented in the JSON (JavaScript Object Notation)
format, which is similar to a Python dictionary.

MISPR uses MongoDB to:

* Add, remove, and search the status of workflows - feature of
  `FireWorks <https://materialsproject.github.io/fireworks/>`__  (required)
* Create computational databases of DFT and MD predicted properties -
  Feature of MISPR (optional but strongly recommended)

Setting up MongoDB
============================
Options for getting MongoDB are:

* Install it yourself locally by following the instructions at
  `MongoDB <https://www.mongodb.com/docs/manual/installation/>`__.
  This is pretty simple and typically works well if you are starting out
  with MISPR and want to learn how to use a database. However, with this
  option, you are limited with the storage space on your local machine and
  you do not have the option to share the database with other users. You
  also need to have the necessary privileges to install mongo on your machine.
* Set up an account using a commercial service, which is typically
  the simplest and easiest to use but is not free of charge for databases
  with large size. Examples of such services include Atlas and MongoDB Atlas,
  which offer 500 MB databases for free. This is typically enough to get
  started for small projects.
* Self-host a MongoDB server or ask your supercomputing center to offer
  MongoDB hosting. This is more complicated than the other options and
  will require continuous maintenance of the server.

After creating a new database, you need to keep record of your credentials.
These will be used later in setting up the configuration files required
by FireWorks.

.. note::
   MongoDB must be accessible from the computers you are using to run
   the workflows.

Testing your MongoDB connection
================================
**Establishing a Connection to MongoDB Using Pymongo:**

You need to import MongoClient from pymongo and then create a new MongoClient instance.
This instance is used to connect to your MongoDB instance:

.. code-block:: python

    from pymongo import MongoClient

    client = MongoClient("mongodb://localhost:27017/")

In this example, we're connecting to a MongoDB instance that runs on the same machine
(localhost) on port 27017, which is the default port for MongoDB.

**Testing the Connection to MongoDB:**

We can check the connection by listing all the databases:

.. code-block:: python

    print(client.list_database_names())

If the connection is successful, this command will return a list of names of the databases that are present in the
MongoDB instance.

Remember, for you to connect to a MongoDB instance, the MongoDB server needs to be installed and running.
If it's not running on localhost:27017, you will need to provide the appropriate connection string.
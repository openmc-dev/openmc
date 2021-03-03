.. _devguide_docker:

======================
Deployment with Docker
======================

OpenMC can be easily deployed using `Docker <https://www.docker.com/>`_ on any
Windows, Mac or Linux system. With Docker running, execute the following
command in the shell to build a `Docker image`_ called ``debian/openmc:latest``:

.. code-block:: sh

    docker build -t debian/openmc:latest https://github.com/openmc-dev/openmc.git#develop

.. note:: This may take 5 -- 10 minutes to run to completion.

This command will execute the instructions in OpenMC's ``Dockerfile`` to
build a Docker image with OpenMC installed. The image includes OpenMC with
MPICH and parallel HDF5 in the ``/opt/openmc`` directory, and
`Miniconda3 <https://conda.io/miniconda.html>`_ with all of the Python
pre-requisites (NumPy, SciPy, Pandas, etc.) installed. The
`NJOY2016 <https://www.njoy21.io/NJOY2016/>`_ codebase is installed in
``/opt/NJOY2016`` to support full functionality and testing of the
``openmc.data`` Python module. The publicly available nuclear data libraries
necessary to run OpenMC's test suite -- including NNDC and WMP cross sections
and ENDF data --  are in the ``/opt/openmc/data directory``, and the
corresponding :envvar:`OPENMC_CROSS_SECTIONS`,
:envvar:`OPENMC_MULTIPOLE_LIBRARY`, and :envvar:`OPENMC_ENDF_DATA`
environment variables are initialized.

After building the Docker image, you can run the following to see the names of
all images on your machine, including ``debian/openmc:latest``:

.. code-block:: sh

    docker image ls

Now you can run the following to create a `Docker container`_ called
``my_openmc`` based on the ``debian/openmc:latest`` image:

.. code-block:: sh

    docker run -it --name=my_openmc debian/openmc:latest

This command will open an interactive shell running from within the
Docker container where you have access to use OpenMC.

.. note:: The ``docker run`` command supports many
          `options <https://docs.docker.com/engine/reference/commandline/run/>`_
          for spawning containers -- including `mounting volumes`_ from the
          host filesystem -- which many users will find useful.

.. _Docker image: https://docs.docker.com/engine/reference/commandline/images/
.. _Docker container: https://www.docker.com/resources/what-container
.. _options: https://docs.docker.com/engine/reference/commandline/run/
.. _mounting volumes: https://docs.docker.com/storage/volumes/

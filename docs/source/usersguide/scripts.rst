.. _usersguide_scripts:

=======================
Executables and Scripts
=======================

.. _scripts_openmc:

----------
``openmc``
----------

Once you have a model built (see :ref:`usersguide_basics`), you can either run
the openmc executable directly from the directory containing your XML input
files, or you can specify as a command-line argument the directory containing
the XML input files.

.. warning::

    OpenMC models should be treated as code, and it is important to be careful with code from untrusted sources.

For example, if your XML input files are in the directory
``/home/username/somemodel/``, one way to run the simulation would be:

.. code-block:: sh

    cd /home/username/somemodel
    openmc

Alternatively, you could run from any directory:

.. code-block:: sh

    openmc /home/username/somemodel

Note that in the latter case, any output files will be placed in the present
working directory which may be different from
``/home/username/somemodel``. ``openmc`` accepts the following command line
flags:

-c, --volume           Run in stochastic volume calculation mode
-e, --event            Run using event-based parallelism
-g, --geometry-debug   Run in geometry debugging mode, where cell overlaps are
                       checked for after each move of a particle
-n, --particles N      Use *N* particles per generation or batch
-p, --plot             Run in plotting mode
-r, --restart file     Restart a previous run from a state point or a particle
                       restart file
-s, --threads N        Run with *N* OpenMP threads
-t, --track            Write tracks for all particles (up to max_tracks)
-q, --verbosity V      Set the output verbosity to *V*
-v, --version          Show version information
-h, --help             Show help message

.. note:: If you're using the Python API, :func:`openmc.run` is equivalent to
          running ``openmc`` from the command line.

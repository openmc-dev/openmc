import os
import logging
from nose.plugins import Plugin

log = logging.getLogger('nose.plugins.nose_mpi')

class NoseMPI(Plugin):

    mpi_np = "0"
    mpi_exec = "/opt/mpich/3.0.4-gnu/bin/mpiexec"

    def options(self, parser, env=os.environ):
        """Define the command line options for plugin."""
        super(NoseMPI, self).options(parser, env)

        parser.add_option(
            "--mpi-np", dest="mpi_np", default=0,
            help="Number of MPI processors to execute OpenMC executable.")

        parser.add_option(
            "--mpi-exec", dest="mpi_exec",
            default="/opt/mpich/3.0.4-gnu/bin/mpiexec",
            help="Absolute path to mpiexec file.")

    def configure(self, options, conf):
        """Configure plugin based on command line options"""
        super(NoseMPI, self).configure(options, conf)

        try:
            mpi_np = int(options.mpi_np)
            self.enabled = True
            NoseMPI.mpi_np = options.mpi_np
        except:
            self.enabled = False
            return

        if not os.path.exists(options.mpi_exec):
            print 'Need to specify valid mpiexec path.'
            exit()

        NoseMPI.mpi_exec = options.mpi_exec

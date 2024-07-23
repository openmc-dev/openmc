import os
import subprocess

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def install(omp=False, mpi=False, phdf5=False, dagmc=False, libmesh=False, ncrystal=False):

    # Build in debug mode by default with support for MCPL
    pip_command = ['pip', 'install', '.[test,vtk]']
    pip_suffix = ['--config-settings=cmake.args=-DCMAKE_BUILD_TYPE=ON;-DOPENMC_USE_MCPL=ON']

    # Turn off OpenMP if specified
    if not omp:
        pip_suffix.append('-DOPENMC_USE_OPENMP=off')

    # Use MPI wrappers when building in parallel
    if mpi:
        pip_suffix.append('-DOPENMC_USE_MPI=on')

    # Tell CMake to prefer parallel HDF5 if specified
    if phdf5:
        if not mpi:
            raise ValueError('Parallel HDF5 must be used in '
                             'conjunction with MPI.')
        pip_suffix.append('-DHDF5_PREFER_PARALLEL=ON')
    else:
        pip_suffix.append('-DHDF5_PREFER_PARALLEL=OFF')

    if dagmc:
        pip_suffix.append('-DOPENMC_USE_DAGMC=ON')
        pip_suffix.append('-DCMAKE_PREFIX_PATH=~/DAGMC')

    if libmesh:
        pip_suffix.append('-DOPENMC_USE_LIBMESH=ON')
        libmesh_path = os.environ.get('HOME') + '/LIBMESH'
        pip_suffix.append('-DCMAKE_PREFIX_PATH=' + libmesh_path)

    if ncrystal:
        pip_suffix.append('-DOPENMC_USE_NCRYSTAL=ON')
        ncrystal_cmake_path = os.environ.get('HOME') + '/ncrystal_inst/lib/cmake'
        pip_suffix.append(f'-DCMAKE_PREFIX_PATH={ncrystal_cmake_path}')

    # Build in coverage mode for coverage testing
    pip_suffix.append('-DOPENMC_ENABLE_COVERAGE=on')

    pip_command.append(';'.join(pip_suffix))
    # Build and install
    print(' '.join(pip_command))
    subprocess.check_call(pip_command)

def main():
    # Convert Travis matrix environment variables into arguments for install()
    omp = (os.environ.get('OMP') == 'y')
    mpi = (os.environ.get('MPI') == 'y')
    phdf5 = (os.environ.get('PHDF5') == 'y')
    dagmc = (os.environ.get('DAGMC') == 'y')
    ncrystal = (os.environ.get('NCRYSTAL') == 'y')
    libmesh = (os.environ.get('LIBMESH') == 'y')

    # Build and install
    install(omp, mpi, phdf5, dagmc, libmesh, ncrystal)

if __name__ == '__main__':
    main()

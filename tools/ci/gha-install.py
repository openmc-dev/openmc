import os
import shutil
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


def install(omp=False, mpi=False, phdf5=False, dagmc=False, libmesh=False):
    # Create build directory and change to it
    shutil.rmtree('build', ignore_errors=True)
    os.mkdir('build')
    os.chdir('build')

    # Build in debug mode by default
    cmake_cmd = ['cmake', '-DCMAKE_BUILD_TYPE=Debug']

    # Turn off OpenMP if specified
    if not omp:
        cmake_cmd.append('-DOPENMC_USE_OPENMP=off')

    # Use MPI wrappers when building in parallel
    if mpi:
        cmake_cmd.append('-DOPENMC_USE_MPI=on')

    # Tell CMake to prefer parallel HDF5 if specified
    if phdf5:
        if not mpi:
            raise ValueError('Parallel HDF5 must be used in '
                             'conjunction with MPI.')
        cmake_cmd.append('-DHDF5_PREFER_PARALLEL=ON')
    else:
        cmake_cmd.append('-DHDF5_PREFER_PARALLEL=OFF')

    if dagmc:
        cmake_cmd.append('-DOPENMC_USE_DAGMC=ON')
        cmake_cmd.append('-DCMAKE_PREFIX_PATH=~/DAGMC')

    if libmesh:
        cmake_cmd.append('-DOPENMC_USE_LIBMESH=ON')
        libmesh_path = os.environ.get('HOME') + '/LIBMESH'
        cmake_cmd.append('-DCMAKE_PREFIX_PATH=' + libmesh_path)

    # Build in coverage mode for coverage testing
    cmake_cmd.append('-DOPENMC_ENABLE_COVERAGE=on')

    # Build and install
    cmake_cmd.append('..')
    print(' '.join(cmake_cmd))
    subprocess.check_call(cmake_cmd)
    subprocess.check_call(['make', '-j4'])
    subprocess.check_call(['sudo', 'make', 'install'])

def main():
    # Convert Travis matrix environment variables into arguments for install()
    omp = (os.environ.get('OMP') == 'y')
    mpi = (os.environ.get('MPI') == 'y')
    phdf5 = (os.environ.get('PHDF5') == 'y')
    dagmc = (os.environ.get('DAGMC') == 'y')
    libmesh = (os.environ.get('LIBMESH') == 'y')

    # Build and install
    install(omp, mpi, phdf5, dagmc, libmesh)

if __name__ == '__main__':
    main()

import os
import shutil
import subprocess

def install(omp=False, mpi=False, phdf5=False, dagmc=False, libmesh=False, ncrystal=False):
    # Create build directory and change to it
    shutil.rmtree('build', ignore_errors=True)
    os.mkdir('build')
    os.chdir('build')

    # Build in debug mode by default with support for MCPL
    cmake_cmd = ['cmake', '-DCMAKE_BUILD_TYPE=Debug', '-DOPENMC_USE_MCPL=on']

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
        cmake_cmd.append('-DOPENMC_USE_UWUW=ON')
        dagmc_path = os.environ.get('HOME') + '/DAGMC'
        cmake_cmd.append('-DCMAKE_PREFIX_PATH=' + dagmc_path)

    if libmesh:
        cmake_cmd.append('-DOPENMC_USE_LIBMESH=ON')
        libmesh_path = os.environ.get('HOME') + '/LIBMESH'
        cmake_cmd.append('-DCMAKE_PREFIX_PATH=' + libmesh_path)

    if ncrystal:
        cmake_cmd.append('-DOPENMC_USE_NCRYSTAL=ON')
        ncrystal_path = os.environ.get('HOME') + '/NCRYSTAL'
        cmake_cmd.append(f'-DCMAKE_PREFIX_PATH={ncrystal_path}')

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
    ncrystal = (os.environ.get('NCRYSTAL') == 'y')
    libmesh = (os.environ.get('LIBMESH') == 'y')

    # Build and install
    install(omp, mpi, phdf5, dagmc, libmesh, ncrystal)

if __name__ == '__main__':
    main()

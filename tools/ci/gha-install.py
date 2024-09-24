import os
import subprocess
import shutil


def install(omp=False, mpi=False, phdf5=False, dagmc=False, libmesh=False, ncrystal=False):
    # List to store the CMake arguments
    cmake_args = ['-DCMAKE_BUILD_TYPE=Debug', '-DOPENMC_USE_MCPL=on']

    # Turn off OpenMP if specified
    if not omp:
        cmake_args.append('-DOPENMC_USE_OPENMP=off')

    # Use MPI wrappers when building in parallel
    if mpi:
        cmake_args.append('-DOPENMC_USE_MPI=on')

    # Tell CMake to prefer parallel HDF5 if specified
    if phdf5:
        if not mpi:
            raise ValueError('Parallel HDF5 must be used in conjunction with MPI.')
        cmake_args.append('-DHDF5_PREFER_PARALLEL=ON')
    else:
        cmake_args.append('-DHDF5_PREFER_PARALLEL=OFF')

    if dagmc:
        cmake_args.append('-DOPENMC_USE_DAGMC=ON')
        cmake_args.append('-DCMAKE_PREFIX_PATH=~/DAGMC')

    if libmesh:
        cmake_args.append('-DOPENMC_USE_LIBMESH=ON')
        libmesh_path = os.environ.get('HOME') + '/LIBMESH'
        cmake_args.append('-DCMAKE_PREFIX_PATH=' + libmesh_path)

    if ncrystal:
        cmake_args.append('-DOPENMC_USE_NCRYSTAL=ON')
        ncrystal_cmake_path = os.environ.get('HOME') + '/ncrystal_inst/lib/cmake'
        cmake_args.append(f'-DCMAKE_PREFIX_PATH={ncrystal_cmake_path}')

    # Build in coverage mode for coverage testing
    cmake_args.append('-DOPENMC_ENABLE_COVERAGE=on')

    # Set environment variable for SKBUILD
    os.environ['SKBUILD_CMAKE_ARGS'] = ';'.join(cmake_args)

    # Run pip to build and install
    subprocess.check_call(['pip', '-v', 'install', '.[test,vtk,ci]'])

    # Using standard CMake method
    # Create build directory and change to it
    shutil.rmtree('build', ignore_errors=True)
    os.mkdir('build')
    os.chdir('build')

    # Add CMake arguments for standard method
    cmake_cmd = ['cmake', '..'] + cmake_args
    print(' '.join(cmake_cmd))

    # Run CMake and build
    subprocess.check_call(cmake_cmd)
    subprocess.check_call(['make', '-j$(nproc)'])
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

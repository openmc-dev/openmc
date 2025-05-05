import os
import sys
import shutil
import subprocess


def install(omp=False, mpi=False, phdf5=False, dagmc=False, libmesh=False):
    # Create build directory and change to it
    shutil.rmtree('build', ignore_errors=True)
    os.mkdir('build')
    os.chdir('build')

    # Build in debug mode by default with support for MCPL
    if sys.platform == 'win32':
        work_dir = os.environ.get('GITHUB_WORKSPACE')
        cmake_cmd = ['cmake', '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_TOOLCHAIN_FILE='+work_dir+'\\vcpkg\\scripts\\buildsystems\\vcpkg.cmake', '-DVCPKG_TARGET_TRIPLET=x64-windows-static']
    else:
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

    # Build in coverage mode for coverage testing
    if sys.platform != 'win32':
        cmake_cmd.append('-DOPENMC_ENABLE_COVERAGE=on')

    # Build and install
    cmake_cmd.append('..')
    print(' '.join(cmake_cmd))
    subprocess.check_call(cmake_cmd)

    if sys.platform == 'win32':
        subprocess.check_call(['cmake', '--build', '.', '--config=Release'])
        subprocess.check_call(['cmake', '--install', '.', '--config=Release'])
    else:
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

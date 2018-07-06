import os
import shutil
import subprocess

MOAB_VERSION= '5.0'

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


def install(omp=False, mpi=False, phdf5=False, cad=False):
    # Create build directory and change to it
    shutil.rmtree('build', ignore_errors=True)
    os.mkdir('build')
    os.chdir('build')

    # Build in debug mode by default
    cmake_cmd = ['cmake', '-Ddebug=on']

    # Turn off OpenMP if specified
    if not omp:
        cmake_cmd.append('-Dopenmp=off')

    # Use MPI wrappers when building in parallel
    if mpi:
        os.environ['FC'] = 'mpifort' if which('mpifort') else 'mpif90'
        os.environ['CC'] = 'mpicc'
        os.environ['CXX'] = 'mpicxx'

    # Tell CMake to prefer parallel HDF5 if specified
    if phdf5:
        if not mpi:
            raise ValueError('Parallel HDF5 must be used in '
                             'conjunction with MPI.')
        cmake_cmd.append('-DHDF5_PREFER_PARALLEL=ON')
    else:
        cmake_cmd.append('-DHDF5_PREFER_PARALLEL=OFF')

    if cad:
        build_dagmc()
        
    # Build and install
    cmake_cmd.append('..')
    print(' '.join(cmake_cmd))
    subprocess.check_call(cmake_cmd)
    subprocess.check_call(['make', '-j'])
    subprocess.check_call(['sudo', 'make', 'install'])

def build_dagmc():

    current_dir = os.getcwd()
    home_dir = os.environ['HOME']
    
    os.chdir(home_dir)
    os.mkdir('MOAB')
    os.chdir('MOAB')
    
    clone_cmd = ['git', 'clone', '-b', 'Version'+MOAB_VERSION, 'https://bitbucket.org/fathomteam/moab']
    subprocess.check_call(clone_cmd)

    os.mkdir('build')
    os.chdir('build')

    cmake_cmd = ['cmake', '../moab', '-DENABLE_HDF5=ON', '-DENABLE_TOOLS=ON', '-DCMAKE_INSTALL_PREFIX='+home_dir]
    subprocess.check_call(cmake_cmd)

    subprocess.check_call(['make', '-j'])

    subprocess.check_call(['make', 'test'])

    subprocess.check_call(['make', 'test'])

def main():
    # Convert Travis matrix environment variables into arguments for install()
    omp = (os.environ.get('OMP') == 'y')
    mpi = (os.environ.get('MPI') == 'y')
    phdf5 = (os.environ.get('PHDF5') == 'y')

    cad = (os.environ.get('CAD') == 'y')

    # Build and install
    install(omp, mpi, phdf5, cad)


if __name__ == '__main__':
    main()

def get_ncrystal_libdir():
    """
    Returns NCrystal library directory or None.
    """
    import shutil
    cmd = shutil.which('ncrystal-config')
    if cmd:
        import subprocess
        res = subprocess.run([cmd,'--show','shlibdir'],capture_output=True)
        if res.returncode == 0:
            return res.stdout.decode('utf8').strip()

def runtime_env():
    """

    Perform any dynamic environment modifications needed for openmc to load
    shared libraries etc.
    """

    #Currently this only concerns NCrystal, which needs the directory with the
    #NCrystal library to be added to the relevant lib path variable:
    import platform
    libpathvar = { 'Linux' : 'LD_LIBRARY_PATH',
                 'Darwin' : 'DYLD_LIBRARY_PATH' }.get(platform.system())
    if not libpathvar:
        return

    #Currently we only need this for NCrystal (in principle it is not needed
    #when NCrystal is from Conda or a system package, but it doesn't actually
    #hurt in that case):
    nc_libdir = get_ncrystal_libdir()
    if not nc_libdir:
        return

    import os
    env = os.environ.copy()
    env[libpathvar] = '%s:%s'%(env.get(libpathvar),nc_libdir)
    return env

def main():
    """
    Entry point wrapper for the binary openmc executable, for when OpenMC is
    installed from a Python wheel.
    """
    import subprocess
    import pathlib
    import sys
    f = pathlib.Path(__file__).parent.joinpath('core','bin','openmc')
    a = sys.argv[:]
    a[0] = f
    env = runtime_env()
    rv = subprocess.run( a, env = env )
    raise SystemExit(rv.returncode)

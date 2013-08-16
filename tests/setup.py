from setuptools import setup
setup(name="nose-mpi",
    entry_points = {
        'nose.plugins':['nose_mpi = nose_mpi:NoseMPI']
        },
    install_requires = ['nose']
)

from skbuild import setup

setup(
    name="openmc",
    version="0.14.4",
    include_package_data=True,
    description="OpenMC Monte Carlo Code",
    author='The OpenMC Development Team',
    license="MIT",
    packages=['openmc'],
    python_requires=">=3.8",
    cmake_source_dir='.',
    zip_safe=False,
    install_requires=[
        "matplotlib",
        "numpy",
        "scipy",
        "ipython",
        "matplotlib",
        "uncertainties",
        "lxml",
        "pandas",
        "h5py"],
)

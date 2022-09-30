from skbuild import setup

setup(
    include_package_data=True,
    description="OpenMC Monte Carlo Code",
    author='The OpenMC Development Team',
    license="MIT",
    packages=['openmc'],
    python_requires=">=3.8",
    cmake_source_dir='src',
)

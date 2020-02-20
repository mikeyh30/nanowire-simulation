from setuptools import setup, find_packages

setup(
    name="nanowire-simulation",
    author="Michael Hynes",
    version="0.1.0",
    packages=find_packages(exclude=[]),
    install_requires=[
        "kwant>=1.4.1",
        "numpy>=1.11.0",
        "scipy>=0.17.0",
        "matplotlib>=1.5.1",
        "sympy>=0.7.6",
        "tinyarray>=1.2",
        "tqdm",
        "pyyaml",
        "pandas"
        ],
    entry_points=
    {
        "console_scripts": [
            "nwsetup = src.setup_sim:main",
            "simulate = src.simulation:main"
        ]
    }
)
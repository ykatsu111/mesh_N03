import sys
from setuptools import setup


setup_args = dict(
    name="mesh_N03",
    version="2.1",
    install_requires=[
        "numpy",
        "pygeoj",
        "tqdm"
    ],
    py_modules=(
        "mesh_N03",
    ),
    url="https://github.com/ykatsu111/mesh_N03.git"
)


if sys.version_info.major == 2:
    raise NotImplementedError(
        "python 2.x is not supported."
    )
else:
    setup(
        **setup_args
    )
#!/usr/bin/env python
import setuptools  # noqa: F401 # lgtm [py/unused-import]
from numpy.distutils.core import setup, Extension


ext = [
    Extension(
        name="kmag",
        sources=["src/KMAG2012.f"],
        # f2py_options=["--quiet"],
        # extra_f77_compile_args=["-w"],
    )
]

setup(ext_modules=ext)

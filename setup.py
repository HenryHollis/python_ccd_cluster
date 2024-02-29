#!usr/bin/env python

import os
import platform
import sys
import glob
import numpy


import glob
import shlex
import shutil
import subprocess
import sysconfig

from contextlib import contextmanager
from pathlib import Path
from select import select
from shutil import which
from time import sleep
from typing import List, Iterable, Iterator, Optional, Tuple, TypeVar, Union

###########################################################################

# Check Python's version info and exit early if it is too old
if sys.version_info < (3, 7):
    print("This module requires Python >= 3.7")
    sys.exit(0)

###########################################################################

from setuptools import setup, Extension

try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:
    bdist_wheel = None

if bdist_wheel is not None:
    class bdist_wheel_abi3(bdist_wheel):
        def get_tag(self):
            python, abi, plat = super().get_tag()
            if python.startswith("cp"):
                # on CPython, our wheels are abi3 and compatible back to 3.5
                return "cp38", "abi3", plat

            return python, abi, plat
else:
    bdist_wheel_abi3 = None

should_build_abi3_wheel = (
    bdist_wheel_abi3 and
    platform.python_implementation() == "CPython" and
    sys.version_info >= (3, 8)
)

# Define the extension
macros = []
if should_build_abi3_wheel:
    macros.append(("Py_LIMITED_API", "0x03090000"))

cmdclass = {}

if should_build_abi3_wheel:
    cmdclass["bdist_wheel"] = bdist_wheel_abi3


louvain_ext = Extension('louvain._c_louvain',
                    sources = glob.glob(os.path.join('src', 'louvain', '*.cpp')),
                    py_limited_api=should_build_abi3_wheel,
                    define_macros=macros,
                    include_dirs=['include_louvain', 'build-deps/install/include', numpy.get_include()],
                    library_dirs=['build-deps/install/lib'],
                    runtime_library_dirs=['build-deps/install/lib'],
                     extra_compile_args=['-std=c++11', '-g'] )

leiden_ext = Extension('leidenalg._c_leiden',
                  sources = glob.glob(os.path.join('src', 'leidenalg', '*.cpp')),
                  py_limited_api=should_build_abi3_wheel,
                  define_macros=macros,
                  libraries = ['libleidenalg', 'igraph'],
                  include_dirs=['include', 'build-deps/install/include', numpy.get_include()],
                  library_dirs=['build-deps/install/lib'],
                  runtime_library_dirs=['build-deps/install/lib'],
                  extra_compile_args=['-std=c++11', '-g'] #g for debugging purposes
        )

setup(
    ext_modules = [ louvain_ext,
                   leiden_ext
    ],
    cmdclass=cmdclass
)

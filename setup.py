#!/usr/bin/env python

import os
import platform
from setuptools import setup
from setuptools.extension import Extension
import sysconfig
import pybind11

# All package metadata is contained in `pyproject.toml`.
# Here, we only run code to set up environment variables for compilation of
# the binary extension and register it with the `ext_modules` argument.

PYBIND11_INCLUDE = os.path.join(pybind11.__path__[0], "include")

# Set up libraries and includes
include_dirs = [
    PYBIND11_INCLUDE,
    './src/BEMEWS/_ext',
    './src/BEMEWS/_ext/mstl',
    './src/BEMEWS/_ext/mstl/math2',
    './src/BEMEWS/_ext/mstl/math2/algebra',
    './src/BEMEWS/_ext/mstl/math2/analysis',
    './src/BEMEWS/_ext/mstl/math2/spline',
    './src/BEMEWS/_ext/mstl/physics'
]

if os.name == 'posix':  # macOS or Linux
    LIBDIR = sysconfig.get_config_var('LIBDIR')
    if platform.system() == 'Darwin':
        # Must have LIBOMP_INCLUDE env variable in macOS
        LIBOMP_INCLUDE = os.environ['LIBOMP_INCLUDE']
        include_dirs = [LIBOMP_INCLUDE] + include_dirs
elif os.name == 'nt':  # Windows
    LIBDIR = sysconfig.get_config_var('LIBDEST')

BEMEWS = Extension('BEMEWS._ext',
    define_macros = [
        ('MAJOR_VERSION', '1'),
        ('MINOR_VERSION', '0')
    ],
    include_dirs = include_dirs,
    #
    # libraries = ['stdc++', 'm', 'gomp', 'python3'],
    #
    library_dirs = [LIBDIR],
    extra_compile_args = [
        '-std=c++17',
        '-fPIC',
        '-nostartfiles'
    ],
    #
    # extra_link_args = ['-shared'],
    #
    sources = [
        './src/BEMEWS/_ext/BEMEWS.cpp',
        './src/BEMEWS/_ext/adiabatic_basis.cpp',
        './src/BEMEWS/_ext/eigenvalues.cpp',
        './src/BEMEWS/_ext/flavour_basis.cpp',
        './src/BEMEWS/_ext/input_class.BEMEWS.cpp',
        './src/BEMEWS/_ext/jacobians.cpp',
        './src/BEMEWS/_ext/mixing_angles.cpp',
        './src/BEMEWS/_ext/output.BEMEWS.cpp',
        './src/BEMEWS/_ext/output_matrix.BEMEWS.cpp',
        './src/BEMEWS/_ext/parameters.cpp',
        './src/BEMEWS/_ext/potentials.cpp',
        './src/BEMEWS/_ext/RK.BEMEWS.cpp',
        './src/BEMEWS/_ext/update.BEMEWS.cpp',
        './src/BEMEWS/_ext/mstl/errors2.cpp',
        './src/BEMEWS/_ext/mstl/messages.cpp',
        './src/BEMEWS/_ext/mstl/miscellaneous functions.cpp',
        './src/BEMEWS/_ext/mstl/stdarg2.cpp',
        './src/BEMEWS/_ext/mstl/math2/algebra/column and row vectors.cpp',
        './src/BEMEWS/_ext/mstl/math2/algebra/linear algebra.cpp',
        './src/BEMEWS/_ext/mstl/math2/algebra/mmatrix.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/algorithm3.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/complex2.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/derivative.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/polynomial.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/roots.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/runge kutta.cpp',
        './src/BEMEWS/_ext/mstl/math2/analysis/special functions.cpp',
        './src/BEMEWS/_ext/mstl/math2/spline/discontinuous.cpp',
        './src/BEMEWS/_ext/mstl/math2/spline/interpolation data.cpp',
        './src/BEMEWS/_ext/mstl/physics/units and constants.cpp'
    ])

setup(ext_modules=[BEMEWS])

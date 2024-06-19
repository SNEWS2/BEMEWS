# to use type "sudo python3 setup.EMEWS.py installon the command line

#!/usr/bin/env python
#
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
from setuptools import setup, find_packages
from setuptools.extension import Extension
import sysconfig
import pybind11

#
# Begin setup
#
setup_keywords = dict()
#

#
# Set other keywords for the setup function.
#
# Use entry_points to let `pip` create executable scripts for each target platform.
# See https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
# setup_keywords['entry_points'] = {'console_scripts': ['to_snowglobes = snewpy.to_snowglobes:generate_time_series', ], },
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['zip_safe'] = False
#setup_keywords['test_suite']='snewpy.test.snewpy_test_suite.snewpy_test_suite'

setup_keywords['extras_require'] = {  # Optional
    'dev': ['pytest'],
    'docs':['numpydoc']
}

LIBOMP_INCLUDE = os.environ['LIBOMP_INCLUDE']
PYBIND11_INCLUDE = os.path.join(pybind11.__path__[0], "include")
if os.name == 'posix':  # macOS or Linux
    LIBDIR = sysconfig.get_config_var('LIBDIR')
elif os.name == 'nt':  # Windows
    LIBDIR = sysconfig.get_config_var('LIBDEST')

EMEWS = Extension('EMEWS',
                    define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
                    include_dirs = [LIBOMP_INCLUDE, PYBIND11_INCLUDE, './src/EMEWS', './src/EMEWS/mstl', './src/EMEWS/mstl/math2', './src/EMEWS/mstl/math2/algebra', './src/EMEWS/mstl/math2/analysis', './src/EMEWS/mstl/math2/spline', './src/EMEWS/mstl/physics'],
                    # libraries = ['stdc++', 'm', 'gomp', 'python3'],
                    library_dirs = [LIBDIR],
                    extra_compile_args = ['-std=c++17', '-fPIC', '-nostartfiles'],
                    # extra_link_args = ['-shared'],
                    sources = ['./src/EMEWS/EMEWS.cpp', './src/EMEWS/adiabatic_basis.cpp', './src/EMEWS/eigenvalues.cpp', './src/EMEWS/flavour_basis.cpp', './src/EMEWS/input_class.EMEWS.cpp', './src/EMEWS/jacobians.cpp', './src/EMEWS/mixing_angles.cpp', './src/EMEWS/output.EMEWS.cpp', './src/EMEWS/output_matrix.EMEWS.cpp', './src/EMEWS/parameters.cpp', './src/EMEWS/potentials.cpp', './src/EMEWS/RK.EMEWS.cpp', './src/EMEWS/update.EMEWS.cpp', './src/EMEWS/mstl/errors2.cpp', './src/EMEWS/mstl/messages.cpp', './src/EMEWS/mstl/miscellaneous functions.cpp', './src/EMEWS/mstl/stdarg2.cpp', './src/EMEWS/mstl/math2/algebra/column and row vectors.cpp', './src/EMEWS/mstl/math2/algebra/linear algebra.cpp', './src/EMEWS/mstl/math2/algebra/mmatrix.cpp', './src/EMEWS/mstl/math2/analysis/algorithm3.cpp', './src/EMEWS/mstl/math2/analysis/complex2.cpp', './src/EMEWS/mstl/math2/analysis/derivative.cpp', './src/EMEWS/mstl/math2/analysis/polynomial.cpp', './src/EMEWS/mstl/math2/analysis/roots.cpp', './src/EMEWS/mstl/math2/analysis/runge kutta.cpp', './src/EMEWS/mstl/math2/analysis/special functions.cpp', './src/EMEWS/mstl/math2/spline/discontinuous.cpp', './src/EMEWS/mstl/math2/spline/interpolation data.cpp', './src/EMEWS/mstl/physics/units and constants.cpp'])

setup_keywords['ext_modules'] = [EMEWS]

#
# Run setup command.
#
setup(**setup_keywords)

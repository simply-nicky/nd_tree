import os
import sys
from setuptools import setup, find_namespace_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
import numpy

__version__ = '1.0.0'

extension_args = {'extra_compile_args': ['-fopenmp', '-std=c++20'],
                  'extra_link_args': ['-lgomp'],
                  'library_dirs': ['/usr/local/lib',
                                   os.path.join(sys.prefix, 'lib')],
                  'include_dirs': [numpy.get_include(),
                                   os.path.join(sys.prefix, 'include'),
                                   os.path.join(os.path.dirname(__file__), 'nd_tree/src')]}

extensions = [Pybind11Extension("nd_tree.src.nd_tree",
                                sources=["nd_tree/src/nd_tree.cpp"],
                                define_macros = [('VERSION_INFO', __version__)],
                                **extension_args),]

with open('README.md', 'r') as readme:
    long_description = readme.read()

setup(name='nd_tree',
      cmdclass={"build_ext": build_ext},
      version=__version__,
      author='Nikolay Ivanov',
      author_email="nikolay.ivanov@desy.de",
      long_description=long_description,
      long_description_content_type='text/markdown',
      packages=find_namespace_packages(),
      install_requires=['numpy',],
      ext_modules=extensions,
      extras_require={"test": "pytest"},
      python_requires='>=3.10')

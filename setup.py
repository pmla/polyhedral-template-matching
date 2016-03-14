#!/usr/bin/env python

# NOTE: This setup.py script requires that Python is compiled with a
# C++ compiler, and that this compiler supports C++11.  This is only
# occationally the case!

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from glob import glob
import os

c_src_files = glob('*.c') + glob('svdpolar/*.c')
cpp_src_files = []
for f in ['main.c', 'selftest.c']:
    c_src_files.remove(f)
extra_compile_args = ['-std=c99']

name = 'ptmmodule'

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(os.path.join(numpy.get_include(), 'numpy'))


setup(name=name,
      version='1.0',
      ext_modules=[Extension(name,
                             c_src_files+cpp_src_files,
                             extra_compile_args=extra_compile_args,
                             )
                  ],
      setup_requires=['numpy'],
      cmdclass={'build_ext':build_ext},
      )

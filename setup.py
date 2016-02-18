#!/usr/bin/env python

raise RuntimeError("THIS setup.py FILE HAS NOT BEEN FINISHED: Please type make to build instead.")

from distutils.core import setup, Extension
from glob import glob

c_src_files = glob('*.c')
cpp_src_files = glob('*.cpp')
extra_compile_args = ['-std=c++11']

name = 'polyhedraltemplatematching'



setup(name=name,
      version='1.0',
      ext_modules=[Extension(name,
                             c_src_files+cpp_src_files,
                             extra_compile_args=extra_compile_args,
                             )
                  ],
      )

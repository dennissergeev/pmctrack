#!/usr/bin/env python
# -*- encoding: utf-8

from pathlib import Path
import subprocess as sb
# import warnings

from distutils.core import setup
from distutils.command.build_py import build_py as _build_py

pwd = Path().absolute()
precc = pwd / 'pmctrack/src/_precc'
fortran_modules = ['kind', 'constants', 'vor_partition', 'cf_synop_check',
                   'synop_check', 'min_z', 'steering_wind', 'link2', 'link',
                   'track_check', 'smth', 'tracking_main']
fortran_modules = [precc / f'{mod}.mod' for mod in fortran_modules]
lib_src = pwd / 'pmctrack' / 'core.so'


class build_py(_build_py):
    # Override the build_py class to
    #  (1) make it also compile the f2py shared object
    #  (2) make python module at . the root module called dynlib
    def run(self):
        sb.call("./compile", shell=True)
        super(build_py, self).run()
        self.copy_file(lib_src,
                       Path(self.build_lib) / 'pmctrack/core.so',
                       preserve_mode=True)
        return

    def finalize_options(self):
        self.set_undefined_options('build', ('build_lib', 'build_lib'))
        _build_py.finalize_options(self)
        return


setup(cmdclass={'build_py': build_py},
      name='pmctrack',
      version='0.0.1',
      description='see readme',
      author='Denis Sergeev',
      author_email='dennis.sergeev@gmail.com',
      packages=['pmctrack'],
      package_dir={'pmctrack': 'pmctrack'},
      # py_modules=['test', ],
      # scripts=['bin/dynlib_init.py', ],
      data_files=[
          ('lib', ['pmctrack/core.so']),
          ('include', fortran_modules),
      ]
      )

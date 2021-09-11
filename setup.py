from setuptools import setup, Extension

"""
A minimalist setup is shown.
"""


setup(name='spkmeansmodule',
      version='1.0',
      description='project CAPI',
      ext_modules=[Extension('spkmeansmodule', sources=['spkmeansmodule.c', 'spkmeans.c'])])
from setuptools import Extension, setup

module = Extension("symnmf_c", sources=['symnmf.c', 'symnmfmodule.c'])
setup(name='symnmf_c',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])
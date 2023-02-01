from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'Bridge',
  ext_modules=[
    Extension('bridge',
              sources=['bridge.pyx'],
              extra_compile_args=['-Ofast'],
              language='c')
    ],
  cmdclass = {'build_ext': build_ext}
)
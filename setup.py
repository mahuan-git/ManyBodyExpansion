from setuptools import find_packages
from setuptools import setup

REQUIRED_PACKAGES = [
    'pyscf',
    'numpy',
    'scipy',
    'openfermion',
]


setup(name='mbe',
      version='0.1',
      description='manybody expansion as a quantum chemistry fragmentation scheme',
      author='Ma Huan',
      author_email='mahuan15@mail.ustc.edu.cn',
      packages=find_packages(),
      install_requires=REQUIRED_PACKAGES,
      platforms=['any'],
     )

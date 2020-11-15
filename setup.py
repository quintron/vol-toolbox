import os.path
import shutil
from setuptools import setup

setup(name='voltoolbox',
      version='1.0.0',
      description='Vanilla vol library',
      author='Quintron',
      packages=['voltoolbox'],
      package_dir={'voltoolbox': 'voltoolbox'},
      package_data={'voltoolbox': ['*.pyd']},
      data_files=[('', ['requirements.txt'])])

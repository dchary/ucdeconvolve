####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from setuptools import setup, find_packages
import pathlib
from pkg_resources import parse_requirements

# Read requirements text file to get required packages
with pathlib.Path('requirements.txt').open() as reqsfile:
    reqs = [str(req) for req in parse_requirements(reqsfile)]

# Setup package with params
setup(name='ucdeconvolve',
      version='0.0.6',
      description='Cell Type Deconvolution For Transcriptomic Data',
      url='https://github.com/dchary/ucdeconvolve',
      author='Daniel Charytonowicz',
      author_email='daniel.charytonowicz@icahn.mssm.edu',
      license='GNU GPLv3',
      install_requires=reqs,
      packages=find_packages(),
      package_data={'ucdeconvolve': ['data/metadata.pkl']},
      include_package_data = True,
      zip_safe=False)
####################################################################################################
# # Copyright (C) 2022-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from setuptools import setup
import pathlib
from pkg_resources import parse_requirements

# Read requirements text file to get required packages
with pathlib.Path('requirements.txt').open() as reqsfile:
    reqs = [str(req) for req in parse_requirements(reqsfile)]

# Setup package with params
setup(name='ucdeconvolve',
      version='0.0.1',
      description='Cell Type Deconvolution For Transcriptomic Data',
      url='https://github.com/dchary/ucdeconvolve',
      author='Daniel Charytonowicz',
      author_email='daniel.charytonowicz@icahn.mssm.edu',
      license='MIT',
      install_requires=reqs,
      packages=['ucdeconvolve'],
      zip_safe=False)
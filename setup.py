"""A setuptools based setup module."""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import re
import ast

_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('app_name/__init__.py', 'rb') as f:
    version = str(ast.literal_eval(_version_re.search(
        f.read().decode('utf-8')).group(1)))

setup(
    name='pysatdata',
    version=version,
    description='Python Space Physics Satellite Data Analysis Toolkit',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/zemarchezi/pySatData',
    author='Jose Paulo Marchezi',
    author_email='jose.marchezi@inpe.br, jpmarchezi@gmail.com',
    license='MIT',
    classifiers=['Development Status :: 2 - Pre-Alpha',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3',
                 ],
    keywords='satellite space data tools',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['requests', 'pyspedas', 'pytplot-mpl-temp>=2.1.17',
                      'cdflib>=0.3.20', 'cdasws>=1.7.24', 'netCDF4',
                      'pywavelets', 'astropy','pyqtgraph', 'loguru',
                      'aacgmv2>=2.6.2', 'matplotlib', 'netCDF4',
                      'pandas', 'pytz', 'scipy', 'tqdm',
                      'xarray','urllib3'],
    python_requires='>=3.6',
    include_package_data=True,
)
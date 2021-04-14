"""A setuptools based setup module."""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open


setup(
    name='pysatdata',
    version='0.0.1',
    description='Python Space Physics Satellite Data Analysis Toolkit',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/spedas/pyspedas',
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
    install_requires=['numpy>=1.20.2', 'requests', 'pytplot>=1.7.24',
                      'cdflib>=0.3.20', 'cdasws>=1.7.24', 'netCDF4',
                      'pywavelets', 'pyqtgraph>=0.11.1', 'astropy'],
    python_requires='>=3.6',
    include_package_data=True,
)
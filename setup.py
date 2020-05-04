from __future__ import absolute_import

from setuptools import setup, find_packages
from codecs import open
import six

with open('README-geodata.rst', encoding='utf-8') as f:
	long_description = f.read()

exec(open('geodata/_version.py').read())

setup(
	name='geodata',
	version=__version__,
	author='Jonas Hoersch (FIAS), Tom Brown (FIAS), Gorm Andresen (Aarhus University). Modifications by Michael Davidson (UCSD).',
	author_email='jonas.hoersch@posteo.de',
	description='Library for fetching and converting weather data to power systems data',
	long_description=long_description,
	license='GPLv3',
	packages=find_packages(exclude=['doc', 'test']),
	include_package_data=True,
	install_requires=['numpy',
					  'scipy',
					  'pandas>=0.22',
					  'bottleneck',
					  'numexpr',
					  'xarray>=0.11.2',
					  'netcdf4',
					  'dask[complete]>=0.18.0',
					  'boto3',
					  'toolz',
					  'pyproj',
					  'requests',
					  'matplotlib',
					  'rasterio',
					  'shapely',
					  'progressbar2',
					  'geopandas'],
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
		'Natural Language :: English',
		'Operating System :: OS Independent',
	])

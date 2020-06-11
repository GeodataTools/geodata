## Copyright 2020 Michael Davidson (UCSD), William Honaker. 

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
GEODATA

Geospatial Data Collection and "Pre-Analysis" Tools
"""

from __future__ import absolute_import

from setuptools import setup, find_packages
from codecs import open
import six
import argparse

requirement_list = ['numpy',
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
				  'geopandas']

parser = argparse.ArgumentParser()
parser.add_argument("command", type=str, help="command")
args = parser.parse_args()
if args.command == 'requirements':
	print ("Requirement list: %s" % requirement_list)
	exit()

with open('README.md', encoding='utf-8') as f:
	long_description = f.read()

exec(open('geodata/_version.py').read())

setup(
	name='geodata',
	version=__version__,
	author='Michael Davidson (UCSD), William Honaker',
	author_email='mrdavidson@ucsd.edu',
	description='Geospatial Data Collection and Pre-Analysis Tools',
	long_description=long_description,
	license='GPLv3',
	packages=find_packages(exclude=['doc', 'test']),
	include_package_data=True,
	install_requires=requirement_list,
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
		'Natural Language :: English',
		'Operating System :: OS Independent',
	])

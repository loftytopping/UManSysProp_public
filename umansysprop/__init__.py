# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright 2014 Dave Jones <dave@waveform.org.uk>.
#
# This file is part of umansysprop.
#
# umansysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# umansysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop. If not, see <http://www.gnu.org/licenses/>.
"""An online system for calculating the properties of individual organic
molecules and ensemble mixtures"""
from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')
import sys

__project__ = 'umansysprop'
__version__ = '0.1'
__keywords__ = ['science', 'organic', 'chemistry']
__author__ = 'Dave Jones, Dave Topping'
__author_email__ = '[Dave Jopnes] dave@waveform.org.uk, [Dave Topping] david.topping@manchester.ac.uk'
__url__ = 'http://umansysprop.readthedocs.org/'
__platforms__ = 'ALL'
__requires__ = ['requests',]
__extra_requires__ = {'server': ['openbabel', 'flask', 'flask-wtf', 'wtforms', 'docutils', 'xlsxwriter'],
    'client': [],
    'doc': ['sphinx'],
    }
if sys.version_info[:2] == (3, 2):
    __extra_requires__['doc'].extend([
    # Particular versions are required for Python 3.2 compatibility.
    # The ordering is reverse because that's what easy_install needs...
    'Jinja<2.7',
    'MarkupSafe<0.16',
    ])
__classifiers__ = [
    'Development Status :: 4 - Beta',
    'Environment :: Web Environment',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    ]
__entry_points__ = {
    'console_scripts': [
    'umansyspropd = umansysprop.server:main',
    ],
    }

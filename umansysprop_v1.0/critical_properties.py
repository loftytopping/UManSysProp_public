# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2016 David Topping.
# All Rights Reserved.
# This file is part of umansysprop.
#
# umansysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# umansysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')


from collections import namedtuple


from . import groups
from . import data


CriticalProperties = namedtuple(
    'CriticalProperties',
    ('temperature', 'pressure', 'volume', 'compressibility'))


def joback_and_reid(compound, boiling_point):
    comp = groups.composition(compound)
    sandb_groups = groups.stein_and_brown(compound)

    n = groups.aggregate_matches(sandb_groups, data.JOBACK_TEMPERATURE)
    temperature = boiling_point / (0.584 + (0.965 * n) - (n ** 2))

    n = groups.aggregate_matches(sandb_groups, data.JOBACK_PRESSURE)
    pressure = (0.113 + (0.0032 * (comp['H'] + comp['n-non-H'])) - n) ** -2

    n = groups.aggregate_matches(sandb_groups, data.JOBACK_VOLUME)
    volume = n + 17.5

    # Factor of 10 added into denominator to account for units of pressure and
    # volume
    compressibility = (pressure * volume) / (8.31451 * 10 * temperature)
    return CriticalProperties(
        temperature=temperature,
        pressure=pressure,
        volume=volume,
        compressibility=compressibility)


def nannoolal(compound, boiling_point):
    comp = groups.composition(compound)
    nannoolal_groups1 = groups.nannoolal_primary(compound)
    nannoolal_groups2 = groups.nannoolal_secondary(compound)
    nannoolal_inter   = groups.nannoolal_interactions(compound)

    n1 = groups.aggregate_matches(
        nannoolal_groups1, data.NANNOOLAL_TEMPERATURE_PRIMARY)
    n2 = groups.aggregate_matches(
        nannoolal_groups2, data.NANNOOLAL_TEMPERATURE_SECONDARY)
    n3 = groups.aggregate_interactions(
        compound, nannoolal_inter, data.NANNOOLAL_TEMPERATURE_INTERACTIONS)
    n = n1 + n2 + n3
    temperature = boiling_point * ((1.0 / (0.9889 + n ** 0.8607)) + 0.699)

    n1 = groups.aggregate_matches(
        nannoolal_groups1, data.NANNOOLAL_PRESSURE_PRIMARY)
    n2 = groups.aggregate_matches(
        nannoolal_groups2, data.NANNOOLAL_PRESSURE_SECONDARY)
    n3 = groups.aggregate_interactions(
        compound, nannoolal_inter, data.NANNOOLAL_PRESSURE_INTERACTIONS)
    num = (comp['mass'] ** -0.14041) / 100.0
    denom = (0.00939 + n1 + n2 + n3) ** 2
    pressure = num / denom

    n1 = groups.aggregate_matches(
        nannoolal_groups1, data.NANNOOLAL_VOLUME_PRIMARY)
    n2 = groups.aggregate_matches(
        nannoolal_groups2, data.NANNOOLAL_VOLUME_SECONDARY)
    n3 = groups.aggregate_interactions(
        compound, nannoolal_inter, data.NANNOOLAL_VOLUME_INTERACTIONS)
    num = n1 + n2 + n3
    denom = comp['n-non-H'] ** -0.2266
    volume = (num / denom) + 86.1539

    # Factor of 10 added into denominator to account for units of pressure and
    # volume
    compressibility = (pressure * volume) / (8.31451 * 10 * temperature)
    return CriticalProperties(
        temperature=temperature,
        pressure=pressure,
        volume=volume,
        compressibility=compressibility)


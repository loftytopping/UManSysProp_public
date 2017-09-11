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


import warnings
from math import sqrt, exp, log, log10, isinf

from . import data
from . import groups
from . import activity_coefficient_models as aiomfac
from .forms import smiles
import numpy as np

class PartitioningIterationLimit(Warning):
    """
    Warning raised when the partitioning loop is terminated due to an excessive
    number of iterations
    """


class PartitioningPrecisionLimit(Warning):
    """
    Warning raised when then partitioning loop is terminated due to the
    precision limit being reached
    """


class IdealPartitioning(object):
    """
    Associates a SMILES compound with various data necessary for the ideal
    partitioning model.
    """

    def __init__(self, compound, abundance, temperature, pressure):
        self.compound = compound
        self.abundance = abundance
        self.molar_mass = compound.molwt
        # Mole-based Ci*, µmol/m³
        self.c_i_star = (1e6 * 10 ** pressure) / (0.000082057 * temperature)
        self.activity_coefficient = 1.0
        # To be calculated iteratively by update()
        self.condensed_abundance = None
        self.activity = None

    def update(self, coa, coefficient):
        self.activity_coefficient=np.clip(coefficient, 1.0e-3, 10000.0)
        self.condensed_abundance = ((1 + self.c_i_star*self.activity_coefficient / coa) ** -1) * self.abundance * (1e12 / 6.023e23)
        self.activity = self.activity_coefficient * (self.condensed_abundance / coa)


class NonIdealPartitioning(object):
    """
    Extends :class:`IdealPartitioning` with the extra data required by the
    AIOMFAC non-ideal partitioning model.
    """

    def __init__(self, compound, abundance, temperature, pressure):
        super(NonIdealPartitioning, self).__init__(
            compound, abundance, temperature, pressure)
        m = groups.aiomfac(compound)
        # XXX What about inorganic?
        #m.update(aiomfac_inorganic(inorganic_ions))

        # XXX What do any of these variables mean?! Unfortunately the former
        # coder left no clue beyond their name ...
        self.q_k_i = { group: count * data.AIOMFAC_QI[group] for group, count in m.items() }
        self.q_i = groups.aggregate_matches(m, data.AIOMFAC_QI)
        self.r_i = groups.aggregate_matches(m, data.AIOMFAC_RI)
        #self.t_i_m_n =

        non_zero_groups = {
            group
            for compound, matches in m.items()
            for group, count in matches.items()
            if count > 0
            }

        q_k_i = {
            compound: {
                group: count * data.AIOMFAC_QI[group]
                for group, count in matches.items()
                }
            for compound, matches in m.items()
            }

        q_i = {
            compound: groups.aggregate_matches(matches, data.AIOMFAC_QI)
            for compound, matches in m.items()
            }

        r_i = {
            compound: groups.aggregate_matches(matches, data.AIOMFAC_RI)
            for compound, matches in m.items()
            }

        t_i_m_n = {
            compound: {
                group1: {
                    group2: q_k_i[compound].get(group1, 0) * exp(
                        -data.AIOMFAC_SR_INTERACTIONS[main_group1][main_group2] /
                        temperature)
                    for group2 in non_zero_groups
                    for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
                    }
                for group1 in non_zero_groups
                for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
                }
            for compound, matches in m.items()
            }

        u_i_m_n = {
            compound: {
                group1: {
                    group2: q_k_i[compound].get(group1, 0) * exp(
                        -data.AIOMFAC_SR_INTERACTIONS[main_group2][main_group1] /
                        temperature)
                    for group2 in non_zero_groups
                    for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
                    }
                for group1 in non_zero_groups
                for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
                }
            for compound, matches in m.items()
            }

        s_i_n = {
            compound: {
                group1: sum(t_i_m_n[compound][group2][group1] for group2 in non_zero_groups)
                for group1 in non_zero_groups
                }
            for compound in t_i_m_n
            }


class WaterMixin(object):
    """
    Water is special-cased in both partition models. This mixin class can be
    added to either of the partitioning classes above to provide water data
    for the model.
    """
    def __init__(self, temperature, humidity):
        # Saturation vapour pressure of water, Pa, Hyland, R. W. and A. Wexler,
        # ASHRAE Trans, 89(2A), 500-519, 1983.
        sat_vap_water = exp(
            (-0.58002206e4 / temperature) + 0.13914993e1 -
            (0.48640239e-1 * temperature) +
            (0.41764768e-4 * temperature ** 2) -
            (0.14452093e-7 * temperature ** 3) +
            (0.65459673e1 * log(temperature)))
        super(WaterMixin, self).__init__(
            compound=smiles('O'),
            abundance=(
                (humidity / 100.0) * 6.023e23 * sat_vap_water * 1e-6 /
                (8.314472 * temperature)
                ),
            temperature=temperature,
            pressure=log10(sat_vap_water * 9.86923267e-6),
            )
        # XXX This is wrong
        self.molar_mass=18.016
        # XXX Not sure why we're overriding the c_i_star calculation here,
        # but it's necessary to match the original...
        self.c_i_star=(1e6 * sat_vap_water) / (8.314472 * temperature)


class WaterIdealPartitioning(WaterMixin, IdealPartitioning):
    pass


class WaterNonIdealPartitioning(WaterMixin, NonIdealPartitioning):
    pass


def aiomfac_organic(compounds):
    return {
        compound: groups.aiomfac(compound)
        for compound in compounds
        }


def aiomfac_inorganic(compounds):
    return {
        compound: groups.matches(data.AIOMFAC_ION_SMARTS, compound)
        for compound in compounds
        }


def aiomfac_salts(compounds):
    m = groups.all_matches(data.AIOMFAC_ION_SMARTS, compounds)
    total_ions = groups.aggregate_matches(m, data.AIOMFAC_ION_CHARGE_ABS)
    cations = {
        ion
        for ion, count in m.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    anions  = {
        ion
        for ion, count in m.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    result = {}
    for cation in cations:
        for anion in anions:
            salt = data.AIOMFAC_ION_SALT[(cation, anion)]
            quantity = 2.0 * m[cation] * m[anion] * sqrt(
                    (data.AIOMFAC_ION_CHARGE_ABS[cation] * data.AIOMFAC_ION_CHARGE_ABS[anion]) /
                    (data.AIOMFAC_SALT_CATION_STOICH[salt] * data.AIOMFAC_SALT_ANION_STOICH[salt])
                    ) / total_ions
            if quantity > 0.0:
                result[salt] = quantity
    return result


def activity_coefficients_sr(organic_compounds, inorganic_ions, temperature):
    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    # XXX What do any of these variables mean?! Unfortunately the former coder
    # left no clue beyond their name ...

    # XXX Need to exclude any groups which have zero matches in the following
    # calculation; the original does this as an "optimization" and yet it
    # produces a different answer. Hmmm ...

    non_zero_groups = {
        group
        for compound, matches in m.items()
        for group, count in matches.items()
        if count > 0
        }

    q_k_i = {
        compound: {
            group: count * data.AIOMFAC_QI[group]
            for group, count in matches.items()
            }
        for compound, matches in m.items()
        }

    q_i = {
        compound: groups.aggregate_matches(matches, data.AIOMFAC_QI)
        for compound, matches in m.items()
        }

    r_i = {
        compound: groups.aggregate_matches(matches, data.AIOMFAC_RI)
        for compound, matches in m.items()
        }

    t_i_m_n = {
        compound: {
            group1: {
                group2: q_k_i[compound].get(group1, 0) * exp(
                    -data.AIOMFAC_SR_INTERACTIONS[main_group1][main_group2] /
                    temperature)
                for group2 in non_zero_groups
                for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
                }
            for group1 in non_zero_groups
            for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
            }
        for compound, matches in m.items()
        }

    u_i_m_n = {
        compound: {
            group1: {
                group2: q_k_i[compound].get(group1, 0) * exp(
                    -data.AIOMFAC_SR_INTERACTIONS[main_group2][main_group1] /
                    temperature)
                for group2 in non_zero_groups
                for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
                }
            for group1 in non_zero_groups
            for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
            }
        for compound, matches in m.items()
        }

    s_i_n = {
        compound: {
            group1: sum(t_i_m_n[compound][group2][group1] for group2 in non_zero_groups)
            for group1 in non_zero_groups
            }
        for compound in t_i_m_n
        }

    water = [c for c in organic_compounds if str(c).strip() == 'O'][0]


def partition_model(organic_compounds, inorganic_ions, vapour_pressures,
        soluble_core, temperature, humidity, ideality):
    water_class = WaterIdealPartitioning #(WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = IdealPartitioning # (NonIdealPartitioning, IdealPartitioning)[ideality]
    m = [water_class(temperature, humidity)]
    m.extend([
        generic_class(
            compound=compound,
            abundance=abundance,
            temperature=temperature,
            pressure=vapour_pressures[compound],
            )
        for compound, abundance in organic_compounds.items()
        ])

    salts = aiomfac_salts(inorganic_ions)
    core_ion = sum(inorganic_ions.values())

    Coa = 1.5
    BCoa = 0.3
    iteration = 1

    #It looks like we need to split the ideal and non-ideal testcases here.
    #1) For ideal runs

    if ideality=='ideal':
        while abs((Coa - BCoa) / Coa) > 1e-8:
            for c in m:
                c.update(Coa, 1.0)
            BCoa = sum(c.condensed_abundance for c in m) + soluble_core.concentration + core_ion
            Coa = (Coa + BCoa * 99.0) / 100.0
            iteration += 1
            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'partition_model iteration limit (200) reached'))
                break
            if isinf(abs((Coa - BCoa) / BCoa)):
                warnings.warn(
                    PartitioningPrecisionLimit(
                        'partition_model precision limit reached'))
                break

        #water, everything_else = m[0], m[1:]

    #2) For nonideal runs

    elif ideality=='nonideal':
        #First get an estimate for the initial activity coefficient calculation
        for c in m:
            c.update(Coa, 1.0)
        organic_compounds_act = {
            c.compound: c.condensed_abundance
            for c in m
            }
        Activity_coefficients_sr=aiomfac.aiomfac_sr(organic_compounds_act, inorganic_ions, temperature)
        Activity_coefficients_mr_lr=aiomfac.aiomfac_mr(organic_compounds_act, inorganic_ions, temperature)
        for c in m:
            c.update(Coa, exp(Activity_coefficients_sr[c.compound]+Activity_coefficients_mr_lr[c.compound]))

        while abs((Coa - BCoa) / Coa) > 1e-8:
            organic_compounds_act = {
                c.compound: c.condensed_abundance
                for c in m
                }
            Activity_coefficients_sr=aiomfac.aiomfac_sr(organic_compounds_act, inorganic_ions, temperature)
            Activity_coefficients_mr_lr=aiomfac.aiomfac_mr(organic_compounds_act, inorganic_ions, temperature)
            for c in m:
                c.update(Coa, exp(Activity_coefficients_sr[c.compound]+Activity_coefficients_mr_lr[c.compound]))
            BCoa = sum(c.condensed_abundance for c in m) + soluble_core.concentration + core_ion
            Coa = (Coa + BCoa * 99.0) / 100.0
            iteration += 1
            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'partition_model iteration limit (200) reached'))
                break
            if isinf(abs((Coa - BCoa) / BCoa)):
                warnings.warn(
                    PartitioningPrecisionLimit(
                        'partition_model precision limit reached'))
                break
        #water, everything_else = m[0], m[1:]
    water, everything_else = m[0], m[1:]

    return water, everything_else

def partition_model_org(organic_compounds, vapour_pressures,
        soluble_core, temperature, humidity, ideality):
    water_class = WaterIdealPartitioning #(WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = IdealPartitioning # (NonIdealPartitioning, IdealPartitioning)[ideality]
    m = [water_class(temperature, humidity)]
    m.extend([
        generic_class(
            compound=compound,
            abundance=abundance,
            temperature=temperature,
            pressure=vapour_pressures[compound],
            )
        for compound, abundance in organic_compounds.items()
        ])

    #salts = aiomfac_salts(inorganic_ions)
    #core_ion = sum(inorganic_ions.values())

    Coa = 1.5
    BCoa = 0.3
    iteration = 1

    #It looks like we need to split the ideal and non-ideal testcases here.
    #1) For ideal runs

    if ideality=='ideal':
        while abs((Coa - BCoa) / Coa) > 1e-8:
            for c in m:
                c.update(Coa, 1.0)
            BCoa = sum(c.condensed_abundance for c in m) + soluble_core.concentration
            Coa = (Coa + BCoa * 99.0) / 100.0
            iteration += 1
            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'partition_model iteration limit (200) reached'))
                break
            if isinf(abs((Coa - BCoa) / BCoa)):
                warnings.warn(
                    PartitioningPrecisionLimit(
                        'partition_model precision limit reached'))
                break

        #water, everything_else = m[0], m[1:]

    #2) For nonideal runs

    elif ideality=='nonideal':
        #First get an estimate for the initial activity coefficient calculation
        for c in m:
            c.update(Coa, 1.0)
        organic_compounds_act = {
            c.compound: c.condensed_abundance
            for c in m
            }
        Activity_coefficients_sr=aiomfac.aiomfac_sr(organic_compounds_act, {}, temperature)
        #Activity_coefficients_mr_lr=aiomfac.aiomfac_mr(organic_compounds_act, inorganic_ions, temperature)
        for c in m:
            c.update(Coa, exp(Activity_coefficients_sr[c.compound]))

        while abs((Coa - BCoa) / Coa) > 1e-8:
            organic_compounds_act = {
                c.compound: c.condensed_abundance
                for c in m
                }
            Activity_coefficients_sr=aiomfac.aiomfac_sr(organic_compounds_act, {}, temperature)
            #Activity_coefficients_mr_lr=aiomfac.aiomfac_mr(organic_compounds_act, inorganic_ions, temperature)
            for c in m:
                c.update(Coa, exp(Activity_coefficients_sr[c.compound]))
            BCoa = sum(c.condensed_abundance for c in m) + soluble_core.concentration
            Coa = (Coa + BCoa * 99.0) / 100.0
            iteration += 1
            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'partition_model iteration limit (200) reached'))
                break
            if isinf(abs((Coa - BCoa) / BCoa)):
                warnings.warn(
                    PartitioningPrecisionLimit(
                        'partition_model precision limit reached'))
                break
        #water, everything_else = m[0], m[1:]
    water, everything_else = m[0], m[1:]

    return water, everything_else



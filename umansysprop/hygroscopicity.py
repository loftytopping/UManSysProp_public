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
        self.activity_coefficient=coefficient
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

def water_activity_difference(organic_compounds, inorganic_ions, ideality, temperature, humidity, Kelvin_factor):
    water = [c for c, value in organic_compounds.items() if str(c).strip() == 'O'][0]
    list_compounds={}
    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    total_abundance = sum(value for key, value in list_compounds.iteritems())
    x_i = {}
    x_i = {
        compound: value/total_abundance for compound, value in list_compounds.items()
        }
    if ideality=='nonideal':
        Activity_coefficients_sr=aiomfac.aiomfac_sr(organic_compounds, inorganic_ions, temperature)
        Activity_coefficients_mr_lr=aiomfac.aiomfac_mr(organic_compounds, inorganic_ions, temperature)
        water_activity=x_i[water]*exp(Activity_coefficients_sr[water]+Activity_coefficients_mr_lr[water])*Kelvin_factor
        org_activity_coeffs = {
            compound :  exp(Activity_coefficients_sr[compound]+Activity_coefficients_mr_lr[compound])
            for compound, value in organic_compounds.items()
            }
    elif ideality=='ideal':
        water_activity=x_i[water]*Kelvin_factor
        org_activity_coeffs = {
            compound :  1.0
            for compound, value in organic_compounds.items()
            }

    difference=water_activity-(humidity/100.0)

    return difference, water_activity, x_i, org_activity_coeffs

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
    salt_density = {}
    salt_mw = {}
    for cation in cations:
        for anion in anions:
            salt = data.AIOMFAC_ION_SALT[(cation, anion)]
            density = data.AIOMFAC_SALT_DENSITY[salt]
            molw = data.AIOMFAC_SALT_MASS[salt]
            quantity = 2.0 * m[cation] * m[anion] * sqrt(
                    (data.AIOMFAC_ION_CHARGE_ABS[cation] * data.AIOMFAC_ION_CHARGE_ABS[anion]) /
                    (data.AIOMFAC_SALT_CATION_STOICH[salt] * data.AIOMFAC_SALT_ANION_STOICH[salt])
                    ) / total_ions
            if quantity > 0.0:
                result[salt] = quantity
                salt_density[salt] = density
                salt_mw[salt]= molw
    return result, salt_density, salt_mw

def growth_factor_model_inorg(inorganic_ions, temperature, humidity, ideality, size_dependance):
    water_class = WaterIdealPartitioning #(WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = IdealPartitioning # (NonIdealPartitioning, IdealPartitioning)[ideality]
    m = [water_class(temperature, humidity)]
    #We no longer need the information below for just inorganics
    #m.extend([
    #    generic_class(
    #        compound=compound,
    #        abundance=abundance,
    #        temperature=temperature,
    #        pressure=vapour_pressures[compound],
    #        )
    #    for compound, abundance in organic_compounds.items()
    #    ])
    salts, salt_density, salt_mw = aiomfac_salts(inorganic_ions) #need this for density information
    core_ion = sum(inorganic_ions.values())

    #Now work out a conversion factor to go from molar amount to wet droplet size
    diameter=size_dependance['diameter']
    surface_tension=size_dependance['surface_tension']
    total_mass_dry=sum(value*salt_mw[salt] for salt, value in salts.items())
    density_dry=1.0/(sum(((value*salt_mw[salt])/total_mass_dry)/salt_density[salt] for salt, value in salts.items()))
    conversion_factor=((4.0/3.0)*3.141592*pow(diameter*0.5e-9,3.0)*density_dry)/total_mass_dry
    #Now work out the equilibrium water concentration using a root finding method
    #First specify the minimum and maximum water mole fractions associated with
    #the mixture based on the requested RH
    mole_frac_min=0.5*(humidity/100.0)
    mole_frac_max=(humidity/100.0)+(1-(humidity/100.0))*0.999999

    #Now work out how many 'soluble' moles there are
    total_moles=core_ion
    #Parameters in the secant method, a and b represent moles of water to update relevant arrays
    a=(total_moles*mole_frac_min)/(1.0-mole_frac_min)
    b=(total_moles*mole_frac_max)/(1.0-mole_frac_max)
    tol=1.0e-5
    #Concentrations of water in both cases
    organic_compounds_act = {
                c.compound: c.condensed_abundance
                for c in m
                }
    water_key = [c for c in organic_compounds_act if str(c).strip() == 'O'][0]
    organic_compounds_act_a = organic_compounds_act.copy()
    organic_compounds_act_b = organic_compounds_act.copy()
    organic_compounds_act_a[water_key]=a
    organic_compounds_act_b[water_key]=b

    total_mass_wet_a=sum(value*salt_mw[salt] for salt, value in salts.items())+ sum((key.molwt*value for key, value in organic_compounds_act_a.items()))
    density_wet_a=1.0/(sum(((value*salt_mw[salt])/total_mass_wet_a)/salt_density[salt] for salt, value in salts.items())+\
         ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
    wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
    Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

    total_mass_wet_b=sum(value*salt_mw[salt] for salt, value in salts.items())+ sum((key.molwt*value for key, value in organic_compounds_act_b.items()))
    density_wet_b=1.0/(sum(((value*salt_mw[salt])/total_mass_wet_b)/salt_density[salt] for salt, value in salts.items())+\
         ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
    wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
    Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

    #Now call the function that calculates the difference between water activity and RH
    f1, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_a)
    f2, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_b)


    f=1.0 #Initial value to initiate loop
    iteration=1
    if abs(f1)<tol:
        x=a
    elif abs(f2)<tol:
        x=b
    elif abs(f1)>tol and abs(f2)>tol:
        while (abs(f)>tol): #think we also need to attach a loop count here just in case
            x=(a+b)/2.0
            organic_compounds_act[water_key]=x
            total_mass_wet=sum(value*salt_mw[salt] for salt, value in salts.items())+ sum((key.molwt*value for key, value in organic_compounds_act.items()))
            density_wet=1.0/(sum(((value*salt_mw[salt])/total_mass_wet)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act[water_key])/total_mass_wet)/1.0e3)
            wet_radius=pow(((total_mass_wet*conversion_factor)/((4.0/3.0)*3.141592*density_wet)), 1.0/3.0)
            Kelvin_factor=exp((2.0*(water_key.molwt*1.0e-3/density_wet)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius))

            f, water_activity, x_i, org_activity_coeffs =water_activity_difference(organic_compounds_act, inorganic_ions,ideality, temperature, humidity, Kelvin_factor)

            old=b
            b=(a+b)/2.0
            x=a
            organic_compounds_act_a[water_key]=x
            total_mass_wet_a=sum(value*salt_mw[salt] for salt, value in salts.items())+ sum((key.molwt*value for key, value in organic_compounds_act_a.items()))
            density_wet_a=1.0/(sum(((value*salt_mw[salt])/total_mass_wet_a)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
            wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
            Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

            fa, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_a)

            x=b
            organic_compounds_act_b[water_key]=x
            total_mass_wet_b=sum(value*salt_mw[salt] for salt, value in salts.items())+ sum((key.molwt*value for key, value in organic_compounds_act_b.items()))
            density_wet_b=1.0/(sum(((value*salt_mw[salt])/total_mass_wet_b)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
            wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
            Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

            fb, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_b)

            iteration+=1

            if fa*fb>0.0:
                b=old
                a=(a+b)/2.0

            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'Water activity root finding algorithm reached iteration limit. This is likely due to non-convergeant \
                         behaviour from activity coefficient predictions. Please try a slightly different solution or review the \
                         physical rationale for your choice of conditions'))
                break

    #Now use the final 'x' value to calculate density and thus growth factor. Note that here our 'dry' mass is only inorganic
    growth_factor=pow((total_mass_wet*density_dry)/(total_mass_dry*density_wet), 1.0/3.0)

    #Now kappa values
    temp=pow(growth_factor,3.0)
    Kappa=1.0-temp+((temp-1.0)*(Kelvin_factor/(humidity/100.0)))

    #Now the equilbrium vapour pressures above the solution for organic compounds
    #org_vp = {
    #    compound: x_i[compound]*pow(10,vapour_pressures[compound])
    #    for compound, value in organic_compounds.items()
    #    }

    mass_frac_sol=total_mass_wet/total_mass_dry

    return growth_factor, Kappa, water_activity, mass_frac_sol


def growth_factor_model_org(organic_compounds, pressure, density, temperature, humidity, ideality, size_dependance):
    water_class = WaterIdealPartitioning #(WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = IdealPartitioning # (NonIdealPartitioning, IdealPartitioning)[ideality]
    m = [water_class(temperature, humidity)]
    m.extend([
        generic_class(
            compound=compound,
            abundance=abundance,
            temperature=temperature,
            pressure=pressure[compound],
            )
        for compound, abundance in organic_compounds.items()
        ])
    #salts, salt_density, salt_mw = aiomfac_salts(inorganic_ions) #need this for density information
    #core_ion = sum(inorganic_ions.values())

    #Now work out a conversion factor to go from molar amount to wet droplet size
    diameter=size_dependance['diameter']
    surface_tension=size_dependance['surface_tension']

    total_mass_dry=sum((key.molwt*value for key, value in organic_compounds.items()))
    density_dry=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items()))

    conversion_factor=((4.0/3.0)*3.141592*pow(diameter*0.5e-9,3.0)*density_dry)/total_mass_dry
    #Now work out the equilibrium water concentration using a root finding method
    #First specify the minimum and maximum water mole fractions associated with
    #the mixture based on the requested RH
    mole_frac_min=0.5*(humidity/100.0)
    mole_frac_max=(humidity/100.0)+(1-(humidity/100.0))*0.999999

    #Now work out how many 'soluble' moles there are
    total_moles=sum((value for key, value in organic_compounds.items()))
    #Parameters in the secant method, a and b represent moles of water to update relevant arrays
    a=(total_moles*mole_frac_min)/(1.0-mole_frac_min)
    b=(total_moles*mole_frac_max)/(1.0-mole_frac_max)
    tol=1.0e-5
    #Concentrations of water in both cases
    organic_compounds_act = {
                c.compound: c.abundance
                for c in m
                }
    water_key = [c for c in organic_compounds_act if str(c).strip() == 'O'][0]
    organic_compounds_act_a = organic_compounds_act.copy()
    organic_compounds_act_b = organic_compounds_act.copy()
    organic_compounds_act_a[water_key]=a
    organic_compounds_act_b[water_key]=b

    total_mass_wet_a=sum((key.molwt*value for key, value in organic_compounds_act_a.items()))
    density_wet_a=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
         ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
    wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
    Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

    total_mass_wet_b=sum((key.molwt*value for key, value in organic_compounds_act_b.items()))
    density_wet_b=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
         ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
    wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
    Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

    #Now call the function that calculates the difference between water activity and RH
    f1, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, {}, ideality, temperature, humidity, Kelvin_factor_a)
    f2, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, {}, ideality, temperature, humidity, Kelvin_factor_b)


    f=1.0 #Initial value to initiate loop
    iteration=1
    if abs(f1)<tol:
        x=a
    elif abs(f2)<tol:
        x=b
    elif abs(f1)>tol and abs(f2)>tol:
        while (abs(f)>tol): #think we also need to attach a loop count here just in case
            x=(a+b)/2.0
            organic_compounds_act[water_key]=x
            total_mass_wet=sum((key.molwt*value for key, value in organic_compounds_act.items()))
            density_wet=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                ((water_key.molwt*organic_compounds_act[water_key])/total_mass_wet)/1.0e3)
            wet_radius=pow(((total_mass_wet*conversion_factor)/((4.0/3.0)*3.141592*density_wet)), 1.0/3.0)
            Kelvin_factor=exp((2.0*(water_key.molwt*1.0e-3/density_wet)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius))

            f, water_activity, x_i, org_activity_coeffs =water_activity_difference(organic_compounds_act, {},ideality, temperature, humidity, Kelvin_factor)

            old=b
            b=(a+b)/2.0
            x=a
            organic_compounds_act_a[water_key]=x
            total_mass_wet_a=sum((key.molwt*value for key, value in organic_compounds_act_a.items()))
            density_wet_a=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
            wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
            Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

            fa, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, {},ideality, temperature, humidity, Kelvin_factor_a)

            x=b
            organic_compounds_act_b[water_key]=x
            total_mass_wet_b=sum((key.molwt*value for key, value in organic_compounds_act_b.items()))
            density_wet_b=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
            wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
            Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

            fb, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, {},ideality, temperature, humidity, Kelvin_factor_b)

            iteration+=1

            if fa*fb>0.0:
                b=old
                a=(a+b)/2.0

            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'Water activity root finding algorithm reached iteration limit. This is likely due to non-convergeant \
                         behaviour from activity coefficient predictions. Please try a slightly different solution or review the \
                         physical rationale for your choice of conditions'))
                break

    #Now use the final 'x' value to calculate density and thus growth factor. Note that here our 'dry' mass is only inorganic
    growth_factor=pow((total_mass_wet*density_dry)/(total_mass_dry*density_wet), 1.0/3.0)

    #Now kappa values
    temp=pow(growth_factor,3.0)
    Kappa=1.0-temp+((temp-1.0)*(Kelvin_factor/(humidity/100.0)))

    #Now the equilbrium vapour pressures above the solution for organic compounds
    org_vp = {
        compound: log10(x_i[compound]*pow(10,pressure[compound])*org_activity_coeffs[compound])
        for compound, value in organic_compounds.items()
        }

    mass_frac_sol=total_mass_wet/total_mass_dry

    return growth_factor, Kappa, water_activity, org_vp, org_activity_coeffs, mass_frac_sol



def growth_factor_model_inorg_org(organic_compounds, inorganic_ions, pressure, density, temperature, humidity, ideality, size_dependance):
    water_class = WaterIdealPartitioning #(WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = IdealPartitioning # (NonIdealPartitioning, IdealPartitioning)[ideality]
    m = [water_class(temperature, humidity)]
    m.extend([
        generic_class(
            compound=compound,
            abundance=abundance,
            temperature=temperature,
            pressure=pressure[compound],
            )
        for compound, abundance in organic_compounds.items()
        ])
    salts, salt_density, salt_mw = aiomfac_salts(inorganic_ions) #need this for density information
    core_ion = sum(inorganic_ions.values())

    #Now work out a conversion factor to go from molar amount to wet droplet size
    diameter=size_dependance['diameter']
    surface_tension=size_dependance['surface_tension']
    #    total_mass_dry=sum(value*salt_mw[salt] for salt, value in salts.items())
    #density_dry=1.0/(sum(((value*salt_mw[salt])/total_mass_dry)/salt_density[salt] for salt, value in salts.items()))

    salt_dens_temp=sum(salt_density[salt] for salt, value in salts.items())
    total_mass_dry=sum((key.molwt*value for key, value in organic_compounds.items()))+sum(value*salt_mw[salt] for salt, value in salts.items())
    density_dry=1.0/(sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
        sum(((value*salt_mw[salt])/total_mass_dry)/salt_density[salt] for salt, value in salts.items()))
        
    density_dry_start=density_dry
    

    conversion_factor=((4.0/3.0)*3.141592*pow(diameter*0.5e-9,3.0)*density_dry)/total_mass_dry
    #Now work out the equilibrium water concentration using a root finding method
    #First specify the minimum and maximum water mole fractions associated with
    #the mixture based on the requested RH
    mole_frac_min=0.5*(humidity/100.0)
    mole_frac_max=(humidity/100.0)+(1-(humidity/100.0))*0.999999

    #Now work out how many 'soluble' moles there are
    total_moles=sum((value for key, value in organic_compounds.items()))+core_ion
    #Parameters in the secant method, a and b represent moles of water to update relevant arrays
    a=(total_moles*mole_frac_min)/(1.0-mole_frac_min)
    b=(total_moles*mole_frac_max)/(1.0-mole_frac_max)
    tol=1.0e-5
    #Concentrations of water in both cases
    organic_compounds_act = {
                c.compound: c.abundance
                for c in m
                }
    water_key = [c for c in organic_compounds_act if str(c).strip() == 'O'][0]
    org_key = [c for c in organic_compounds_act if str(c).strip() != 'O'][0]
    organic_compounds_act_a = organic_compounds_act.copy()
    organic_compounds_act_b = organic_compounds_act.copy()
    organic_compounds_act_a[water_key]=a
    organic_compounds_act_b[water_key]=b
    #There are some values we dont need to calculate twice
    density_inorg_sum=sum(((value*salt_mw[salt])/total_mass_dry)/salt_density[salt] for salt, value in salts.items())
    density_org_sum=sum(((value*key.molwt)/total_mass_dry)/(density[key]*1.0e3) for key, value in organic_compounds.items())

    total_mass_wet_a=sum((key.molwt*value for key, value in organic_compounds_act_a.items()))+\
        sum(value*salt_mw[salt] for salt, value in salts.items())
    density_wet_a=1.0/(sum(((value*key.molwt)/total_mass_wet_a)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
         sum(((value*salt_mw[salt])/total_mass_wet_a)/salt_density[salt] for salt, value in salts.items())+\
        ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
    wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
    Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

    total_mass_wet_b=sum((key.molwt*value for key, value in organic_compounds_act_b.items()))+sum(value*salt_mw[salt] for salt, value in salts.items())
    density_wet_b=1.0/(sum(((value*key.molwt)/total_mass_wet_b)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
         sum(((value*salt_mw[salt])/total_mass_wet_b)/salt_density[salt] for salt, value in salts.items())+\
         ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
    wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
    Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

    #Now call the function that calculates the difference between water activity and RH
    f1, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, inorganic_ions, ideality, temperature, humidity, Kelvin_factor_a)
    f2, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, inorganic_ions, ideality, temperature, humidity, Kelvin_factor_b)

    #NOTE D TOPPING 1/6/2015: CURRENTLY DANGER OF NOT BOUNDING ROOT SEARCH, DEPENDING ON COMPOSITION AND PERFORMANCE OF ACTIVITY MODEL
    #THIS NEEDS CHECKING AND A SAFETY MECHANISM PUT IN PLACE

    f=1.0 #Initial value to initiate loop
    iteration=1
    if abs(f1)<tol:
        x=a
    elif abs(f2)<tol:
        x=b
    elif abs(f1)>tol and abs(f2)>tol:
        while (abs(f)>tol): #think we also need to attach a loop count here just in case
            x=(a+b)/2.0
            organic_compounds_act[water_key]=x
            total_mass_wet=sum((key.molwt*value for key, value in organic_compounds_act.items()))+sum(value*salt_mw[salt] for salt, value in salts.items())
            density_wet=1.0/(sum(((value*key.molwt)/total_mass_wet)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                sum(((value*salt_mw[salt])/total_mass_wet)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act[water_key])/total_mass_wet)/1.0e3)
            wet_radius=pow(((total_mass_wet*conversion_factor)/((4.0/3.0)*3.141592*density_wet)), 1.0/3.0)
            Kelvin_factor=exp((2.0*(water_key.molwt*1.0e-3/density_wet)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius))

            f, water_activity, x_i, org_activity_coeffs =water_activity_difference(organic_compounds_act, inorganic_ions,ideality, temperature, humidity, Kelvin_factor)

            old=b
            b=(a+b)/2.0
            x=a
            organic_compounds_act_a[water_key]=x
            total_mass_wet_a=sum((key.molwt*value for key, value in organic_compounds_act_a.items()))+sum(value*salt_mw[salt] for salt, value in salts.items())
            density_wet_a=1.0/(sum(((value*key.molwt)/total_mass_wet_a)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                sum(((value*salt_mw[salt])/total_mass_wet_a)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act_a[water_key])/total_mass_wet_a)/1.0e3)
            wet_radius_a=pow(((total_mass_wet_a*conversion_factor)/((4.0/3.0)*3.141592*density_wet_a)), 1.0/3.0)
            Kelvin_factor_a=exp((2.0*(water_key.molwt*1.0e-3/density_wet_a)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_a))

            fa, water_activity_a, x_ia, org_activity_coeffs_a = water_activity_difference(organic_compounds_act_a, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_a)

            x=b
            organic_compounds_act_b[water_key]=x
            total_mass_wet_b=sum((key.molwt*value for key, value in organic_compounds_act_b.items()))+sum(value*salt_mw[salt] for salt, value in salts.items())
            density_wet_b=1.0/(sum(((value*key.molwt)/total_mass_wet_b)/(density[key]*1.0e3) for key, value in organic_compounds.items())+\
                sum(((value*salt_mw[salt])/total_mass_wet_b)/salt_density[salt] for salt, value in salts.items())+\
                ((water_key.molwt*organic_compounds_act_b[water_key])/total_mass_wet_b)/1.0e3)
            wet_radius_b=pow(((total_mass_wet_b*conversion_factor)/((4.0/3.0)*3.141592*density_wet_b)), 1.0/3.0)
            Kelvin_factor_b=exp((2.0*(water_key.molwt*1.0e-3/density_wet_b)*surface_tension*1.0e-3)/(8.3142*temperature*wet_radius_b))

            fb, water_activity_b, x_ib, org_activity_coeffs_b = water_activity_difference(organic_compounds_act_b, inorganic_ions,ideality, temperature, humidity, Kelvin_factor_b)

            iteration+=1

            if fa*fb>0.0:
                b=old
                a=(a+b)/2.0

            if iteration > 200:
                warnings.warn(
                    PartitioningIterationLimit(
                        'Water activity root finding algorithm reached iteration limit. This is likely due to non-convergeant \
                         behaviour from activity coefficient predictions. Please try a slightly different solution or review the \
                         physical rationale for your choice of conditions'))
                break

    #Now use the final 'x' value to calculate density and thus growth factor. Note that here our 'dry' mass is only inorganic
    growth_factor=pow((total_mass_wet*density_dry)/(total_mass_dry*density_wet), 1.0/3.0)

    #Now kappa values
    temp=pow(growth_factor,3.0)
    Kappa=1.0-temp+((temp-1.0)*(Kelvin_factor/(humidity/100.0)))

    #Now the equilbrium vapour pressures above the solution for organic compounds
    org_vp = {
        compound: log10(x_i[compound]*pow(10,pressure[compound])*org_activity_coeffs[compound])
        for compound, value in organic_compounds.items()
        }

    mass_frac_sol=total_mass_wet/total_mass_dry

    return growth_factor, Kappa, water_activity, org_vp, org_activity_coeffs, mass_frac_sol
    #return water_activity_b, Kappa, water_activity, org_vp, org_activity_coeffs, mass_frac_sol


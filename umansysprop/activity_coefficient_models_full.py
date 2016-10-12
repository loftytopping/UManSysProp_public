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
from .forms import smiles


#class PartitioningIterationLimit(Warning):
#    """
#    Warning raised when the partitioning loop is terminated due to an excessive
#    number of iterations
#    """


#class PartitioningPrecisionLimit(Warning):
#    """
#    Warning raised when then partitioning loop is terminated due to the
#    precision limit being reached
#    """


class AIOMFAC_class(object):
    """
    Associates a SMILES compound with various data necessary for the ideal
    partitioning model.
    """

    def __init__(self, compound, abundance):
        self.compound = compound
        self.abundance = abundance
        self.molar_mass = compound.molwt

    def update(self, activity_coefficient):
        self.activity_coefficient = activity_coefficient


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

def aiomfac_mr(organic_compounds, inorganic_ions, temperature):
    m = aiomfac_organic(organic_compounds)
    i = aiomfac_inorganic(inorganic_ions)
    full = aiomfac_organic(organic_compounds)
    full.update(aiomfac_inorganic(inorganic_ions))

    #Calculate the solvent molecular weight from all functional groups
    #old code
    #solvent_kg=0
    #for i in range(len(org_amount)):
    #    solvent_kg=solvent_kg+org_amount[i,0]*M_sol[0,i]
    solvent_kg=sum(key[0].molwt*value for key, value in organic_compounds.iteritems())

    #now cycle through all organic functional groups and calculate the mole fraction of each
    #for k in range(max_group_num_org_main):
    #    func_groups_amount[0,k]=numpy.sum(numpy.dot(org_group_stoich_new_main[:,k],org_amount))
    #mole_frac_func_group[0,:]=func_groups_amount[0,:]/(numpy.sum(func_groups_amount))
    #lets follow the mantra of the SR method given below
    total_non_zero_groups=0
    non_zero_groups = {
        group
        for compound, matches in m.items()
        for group, count in matches.items()
        if count > 0
        }

    for compound, matches in m.items():
        for group, count in matches.items():
            if count > 0:
                total_non_zero_groups=total_non_zero_groups+count

    conc_non_zero_groups = {
        groups: sum(m[compound].get(groups, 0) for compound, matches in m.items())
        for groups in non_zero_groups
                           }

    mole_frac_non_zero_groups = {
        groups: conc_non_zero_groups[groups]/total_non_zero_groups
        for groups in non_zero_groups
                                }

    #Now calculate the mole fraction of each organic compound with each other
    #for k in range(org_molecules):
    #    mole_frac_solv[0,k]=org_amount[k,0]/(numpy.sum(org_amount))
    total_abundance_org = sum(value for key, value in organic_compounds.iteritems())
    mole_frac_solv = {
        compound: value/total_abundance_org for compound, value in organic_compounds.iteritems()
        }

    #now calculate properties associated with the ions - molalities
    

def aiomfac_sr(organic_compounds, inorganic_ions, temperature):
    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
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
            #for group, count in non_zero_groups
            }
        for compound, matches in m.items()
        }
    
    #q_k_i = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count

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
    
    #water = [c for c in organic_compounds if str(c).strip() == 'O'][0]

    #Now adding the rest of the calculation to complete predictions of activity 
    #coefficients from the UNIFAC portion of AIOMFAC [sr]
    
    ##mole fraction
    #list_compounds=organic_compounds.update(inorganic_ions)

    
    total_abundance = sum(value for key, value in list_compounds.iteritems())
    x_i = {
        compound: value/total_abundance for compound, value in list_compounds.iteritems()
        }

    #Combinatorial activity coefficient
    Q = sum(x_i[compound]*q_i[compound] for compound, matches in m.items())

    omega_i = {
        compound: log(q_i[compound]/Q) for compound, matches  in m.items()
        }

    R = sum(x_i[compound]*r_i[compound] for compound, matches  in m.items())

    row_i = {
        compound: r_i[compound]/R for compound, matches  in m.items()
        }

    P_i = {
        compound: log(r_i[compound]/R) for compound, matches  in m.items()
        }

    delta_i = {
        compound: row_i[compound]*(5.0*Q-1.0) for compound, matches  in m.items()
        }

    cross_i = {
        compound: 5.0*q_i[compound] for compound, matches  in m.items()
        }

    Ln_gamma_i_C = {
        compound: 1.0+delta_i[compound]+P_i[compound]+cross_i[compound]*(omega_i[compound]-P_i[compound]-1.0) for compound, matches  in m.items()
        }

    #Residual activity coefficient (Ln_gamma_R)
    uu_m_n = {
        group1: {
             group2: sum(x_i[compound]*u_i_m_n[compound][group1][group2] for compound, matches  in m.items())
             for group2 in non_zero_groups
             }
        for group1 in non_zero_groups
        }

    ss_n = {
        group1: sum(x_i[compound]*s_i_n[compound][group1] for compound, matches in m.items())
        for group1 in non_zero_groups
        }

    xx_weird_i_n = {
        compound: {
            group1: log(s_i_n[compound][group1]/ss_n[group1]) 
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }

    brac1 = {
        compound: {
            group1: sum(u_i_m_n[compound][group2][group1]/s_i_n[compound][group2] for group2 in non_zero_groups)      
            for group1 in non_zero_groups
            }
       for compound, matches in m.items()
       }
   
    brac2 = {
        compound: {
            group1: sum(uu_m_n[group2][group1]/ss_n[group2] for group2 in non_zero_groups)
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }
    
    summation2 = {
        compound: {
            group1: xx_weird_i_n[compound][group1]-omega_i[compound]+brac1[compound][group1]-brac2[compound][group1]
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }
            
    #compared to the original code we have taken out the use of 'summation1'
    Ln_gamma_i_R = {
        compound: sum(q_k_i[compound].get(group1, 0)*summation2[compound][group1] for group1 in non_zero_groups)         
        for compound, matches in m.items()
        }       

    #Final activity coefficient    
    Ln_gamma_i_SR = {
        compound:Ln_gamma_i_C[compound]+Ln_gamma_i_R[compound]
        for compound, matches  in m.items()
        }   

    return Ln_gamma_i_SR            
           

def calculate_activities(organic_compounds, inorganic_ions, temperature):
    #water_class = (WaterNonIdealPartitioning, WaterIdealPartitioning)[ideality]
    generic_class = (AIOMFAC_class)
    #m = [water_class(temperature, humidity)]
    #m.extend([
    #    generic_class(
    #        compound=compound,
    #        abundance=abundance,
    #        temperature=temperature,
    #        )
    #    for compound, abundance in organic_compounds.items()
    #    ])
    m=[generic_class(
         compound=compound,
         abundance=abundance
         )
      for compound, abundance in organic_compounds.items()
      ]
    m.extend([generic_class(
         compound=compound,
         abundance=abundance
         )
      for compound, abundance in inorganic_ions.items()
      ])

    #BCoa = sum(c.abundance for c in m)
    salts = aiomfac_salts(inorganic_ions)
    core_ion = sum(inorganic_ions.values())

    #Now calculate activity coefficients associated with each component
    #return it as a value associated with each 

    Activity_coefficients=aiomfac_sr(organic_compounds, inorganic_ions, temperature)
    step=0
    for c in m:
        x=m[step].compound
        c.update(Activity_coefficients[x])
        step+=1

    #for compound, abundance in m.items(0

    return Activity_coefficients, m



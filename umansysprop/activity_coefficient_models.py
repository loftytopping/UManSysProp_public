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

import pdb

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


def aiomfac_sr(organic_compounds, inorganic_ions, temperature):
    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    m_org = aiomfac_organic(organic_compounds)
    m_inorg = aiomfac_inorganic(inorganic_ions)

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

    #Now work out a conversion factor if ions are present
    solvent_kg=sum((key.molwt*1e-3)*value for key, value in organic_compounds.items())
    m_ions = groups.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
    ion_molalities = {
        ions: value/solvent_kg
        for ions, value in m_ions.iteritems()
                     }
    #Cation/Anion molality
    cation_molality = {
        ion: m_ions[ion]/solvent_kg
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    anion_molality = {
        ion: m_ions[ion]/solvent_kg
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    #Ionic strength
    Ionic_strength = 0
    for c, value in cation_molality.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[c], 2)))
    for a, value in anion_molality.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[a], 2)))

    Ionic_strength = Ionic_strength*0.5
    ionic_conversion_factor=0.0
    if (Ionic_strength > 0.0):
            total_abundance_org = sum(value for key, value in organic_compounds.iteritems())
            x_i_org = {
                compound: value/total_abundance_org for compound, value in organic_compounds.items()
                }
            mean_mw_solvent=sum(x_i_org[compound]*(compound.molwt*1e-3) for compound, value in organic_compounds.items())
            sum_molalities=sum(value for ion, value in ion_molalities.items())
            ionic_conversion_factor=log((0.01801528/mean_mw_solvent)+0.01801528*sum_molalities)


    #Final activity coefficient
    Ln_gamma_i_SR = {
        compound:Ln_gamma_i_C[compound]+Ln_gamma_i_R[compound]
        for compound, matches  in m.items()
        }
    #convert scale from mole fraction to molality for ions
    for compound, matches in inorganic_ions.iteritems():
        Ln_gamma_i_SR[compound] = Ln_gamma_i_SR[compound] - ionic_conversion_factor
                

    return Ln_gamma_i_SR

def aiomfac_mr(organic_compounds, inorganic_ions, temperature):

    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    m_org = aiomfac_organic(organic_compounds)
    m_inorg = aiomfac_inorganic(inorganic_ions)

    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    #Mole fraction of all compounds

    non_zero_groups = {
        group
        for compound, matches in m_org.items()
        for group, count in matches.items()
        if count > 0
        }

    ##-------------TESTING AIOMFAC MR TERMS------------------

    solvent_kg=sum((key.molwt*1e-3)*value for key, value in organic_compounds.items())

    total_non_zero_groups=0
    #conc_non_zero_groups={}
    for compound, matches in m_org.items():
        for group, count in matches.items():
            if count > 0:
                #conc_non_zero_groups[group]=conc_non_zero_groups[group]+count
                total_non_zero_groups=total_non_zero_groups+count

    #only over the organic functional groups
    conc_non_zero_groups = {
        groups: sum(m_org[compound].get(groups, 0) for compound, matches in m_org.items())
        for groups in non_zero_groups
                           }
    #again, only over the organic compounds
    mole_frac_non_zero_groups = {
        groups: conc_non_zero_groups[groups]/total_non_zero_groups
        for groups in non_zero_groups
                                }
    #now calculate properties associated with the ions - molalities
    m_ions = groups.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
    ion_molalities = {
        ions: value/solvent_kg
        for ions, value in m_ions.iteritems()
                     }

    #Cation/Anion molality
    cations = {
        ion
        for ion, count in m_ions.items() #ion is a numeric key, not a pybel object
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    cation_molality = {
        ion: m_ions[ion]/solvent_kg
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    anions  = {
        ion
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    anion_molality = {
        ion: m_ions[ion]/solvent_kg
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }

    #Ionic strength
    Ionic_strength = 0
    for c, value in cation_molality.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[c], 2)))
    for a, value in anion_molality.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[a], 2)))

    Ionic_strength = Ionic_strength*0.5


    ##Now generate the interaction parameters needed between ions and organic groups
    ##b1_ki_temp, b2_ki_temp, b3_ki_temp
    #remeber that interaction parameters are given between 'main' organic groups hence
    #need to convert from non_zero_groups to main as in SR

    #First initialise all variables
    b1_ki_temp = {
        groups: {
            ion: 0.0
            for ion, count in m_ions.items() #remember this 'ion' is numeric not a pybel object
            }
        for groups in non_zero_groups
        }

    b2_ki_temp = {
        groups: {
            ion: 0.0
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups
        }

    b3_ki_temp = {
        groups: {
            ion: 1.2 #b3_ki_full *NOTE:assumed to be constant in Zuend et al paper
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups
        }
    b1_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b1[main_group][ion-1]
            for ion, count in m_ions.items() #remember this 'ion' is numeric not a pybel object
            }
        for groups in non_zero_groups
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }

    b2_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b2[main_group][ion-1]
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }
    ##b1ca, b2ca, b3ca
    b1_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }
    b2_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }
    b3_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }
    b1_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b1[cations-1][anions-1]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }
    b2_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b2[cations-1][anions-1]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }
    b3_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b3[cations-1][anions-1]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    ##c1ca, c2ca, c3ca
    c1_ca_temp = {
        cations: {
            anions: 0.0 #data.AIOMFAC_MR_ION_ION_INTERACTIONS_c1[cations][anions]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    c2_ca_temp = {
        cations: {
            anions: 0.0  #populate default value and change in following if needs be
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    c2_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c2[cations-1][anions-1]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    c3_ca_temp = {
        cations: {
            anions: 0.6 #populate default value and change in following if needs be
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    c3_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c3[cations-1][anions-1]
            for anions, value in anion_molality.items()
            }
        for cations, value in cation_molality.items()
        }

    #Cation-cation interaction parameter matrix
    #for clarity, update all possible ions to have 0.0 value so this can be changed
    cation_dict={202:0.0,203:0.0,204:0.0,205:0.0,221:0.0,223:0.0}
    anion_dict={242:0.0,245:0.0,261:0.0}

    Rc_cdash = {}
    Rc_cdash = {
        cations: {
            cations_dash:value
            for cations_dash, value in cation_dict.items()
            }
        for cations, value in cation_dict.items()
        }
    Rc_cdash[204][205]=-0.154486414200000e+00
    Rc_cdash[205][204]=-0.154486414200000e+00
    #Molecular weight of seperate functional groups
    Mk_sol = {
        groups: data.AIOMFAC_MASS[groups]
        for groups in non_zero_groups
        }

    Qc_cdash = {
        cations: {
            cations_dash:{
                anions: 0.0
                for anions, value in anion_dict.items()
                }
            for cations_dash, value in cation_dict.items()
            }
        for cations, value in cation_dict.items()
        }
    Qc_cdash[204][205][248]=0.448354085300000e-03
    Qc_cdash[205][204][248]=0.448354085300000e-03
    #---------Now perform the interaction parameter calculations----------------

    #Organic-ion interactions
    if (Ionic_strength<250.0):
        Bki = {
            groups: {
                ions : b1_ki_temp[groups][ions]+b2_ki_temp[groups][ions]*exp(-1.0*b3_ki_temp[groups][ions]*(pow(Ionic_strength,0.5)))
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups
            }
        Bki_dash = {
            groups: {
                #Bki_dash[k,i]=-0.5*b2_ki_temp[k,i]*b3_ki_temp[k,i]*(numpy.power(I,-0.5))*numpy.exp(-1*b3_ki_temp[k,i]*numpy.power(I,0.5))
                ions: -0.5*b2_ki_temp[groups][ions]*b3_ki_temp[groups][ions]*(pow(Ionic_strength,-0.5))*(exp(-1*b3_ki_temp[groups][ions]*pow(Ionic_strength,0.5)))
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups
            }

        Bca = {
            cations: {
                anions : b1_ca_temp[cations][anions]+b2_ca_temp[cations][anions]*exp(-1.0*b3_ca_temp[cations][anions]*(pow(Ionic_strength,0.5)))
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }
        Bca_dash = {
            cations: {
                #Bca_dash[c,a]=-0.5*b2_ca_temp[c,a]*b3_ca_temp[c,a]*(numpy.power(I,-0.5))*numpy.exp(-1*b3_ca_temp[c,a]*numpy.power(I,0.5))
                anions: -0.5*b2_ca_temp[cations][anions]*b3_ca_temp[cations][anions]*(pow(Ionic_strength,-0.5))*(exp(-1*b3_ca_temp[cations][anions]*pow(Ionic_strength,0.5)))
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }

        Cca = {
            cations: {
                anions : c1_ca_temp[cations][anions]+c2_ca_temp[cations][anions]*exp(-1.0*c3_ca_temp[cations][anions]*(pow(Ionic_strength,0.5)))
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }
        Cca_dash = {
            cations: {
                #Cca_dash[c,a]=-0.5*c2_ca_temp[c,a]*c3_ca_temp[c,a]*(numpy.power(I,-0.5))*numpy.exp(-1*c3_ca_temp[c,a]*numpy.power(I,0.5))
                anions: -0.5*c2_ca_temp[cations][anions]*c3_ca_temp[cations][anions]*(pow(Ionic_strength,-0.5))*(exp(-1*c3_ca_temp[cations][anions]*pow(Ionic_strength,0.5)))
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }

    elif (Ionic_strength>=250.0):
        Bki = {
            groups: {
                ions : b1_ki_temp[groups][ions]
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups
            }
        Bki_dash = {
            groups: {
                ions: 0.0
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups
            }

        Bca = {
            cations: {
                anions : b1_ca_temp[cations][anions]
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }
        Bca_dash = {
            cations: {
                anions: 0.0
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }

        Cca = {
            cations: {
                anions : c1_ca_temp[cations][anions]
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }
        Cca_dash = {
            cations: {
                anions: 0.0
                for anions, value in anion_molality.items()
                }
            for cations, value in cation_molality.items()
            }

    #Average molecular weight of the functional groups
    M_solv_mix=sum(Mk_sol[key]*mole_frac_non_zero_groups[key] for key in non_zero_groups)
    #Now calculate the activity coefficient of the separate organic functional groups
    summation1 = {
        groups: sum (Bki[groups][ions]*value for ions, value in ion_molalities.iteritems())
        for groups in non_zero_groups
        }
    temp_summation1=sum((Bki[groups][ions]+Ionic_strength*Bki_dash[groups][ions])*\
            (mole_frac_non_zero_groups[groups]*value)
            for ions, value in ion_molalities.iteritems()
            for groups in non_zero_groups)

    summation2 = {
        groups:temp_summation1*(Mk_sol[groups]/M_solv_mix)
        for groups in non_zero_groups
        }
    #A2) ion-ion interactions
    temp_summation2=sum((Bca[cations][anions]+Ionic_strength*Bca_dash[cations][anions])*cation_molality[cations]*anion_molality[anions]
        for anions, value in anion_molality.items()
        for cations, value in cation_molality.items())
    summation3 = {
        groups:temp_summation2*Mk_sol[groups]
        for groups in non_zero_groups
        }
    #A3) ion-ion interactions
    temp_summation3 = 0
    for c, value in cation_molality.items():
        temp_summation3 = temp_summation3 + (value*(abs(data.AIOMFAC_ION_CHARGE[c])))
    for a, value in anion_molality.items():
        temp_summation3 = temp_summation3 + (value*(abs(data.AIOMFAC_ION_CHARGE[a])))
    temp_summation4 = sum ((2.0*Cca[cations][anions]+Ionic_strength*Cca_dash[cations][anions])*ion_molalities[cations]*ion_molalities[anions]
        for anions, value in anion_molality.items()
        for cations, value in cation_molality.items())
    summation4 = {
        groups:temp_summation3*temp_summation4*Mk_sol[groups]
        for groups in non_zero_groups
        }
    #A4) cation, cation (only need to do if more than one cation)
    temp_summation5 = sum (Rc_cdash[cations][cations_dash]*ion_molalities[cations]*ion_molalities[cations_dash]
        for cations_dash, value in cation_molality.items()
        for cations, value in cation_molality.items())
    summation5 = {
        groups:temp_summation5*Mk_sol[groups]
        for groups in non_zero_groups
        }
    temp_summation6 = sum(2.0*Qc_cdash[cations][cations_dash][anions]*c_value*cdash_value*a_value
        for anions, a_value in anion_molality.items()
        for cations_dash, cdash_value in cation_molality.items()
        for cations, c_value in cation_molality.items())
    summation6 = {
        groups:temp_summation6*Mk_sol[groups]
        for groups in non_zero_groups
        }
    ln_gamma_k_MR = {
        groups: summation1[groups]-summation2[groups]-summation3[groups]-summation4[groups]-summation5[groups]-summation6[groups]
        for groups in non_zero_groups
        }

    #Ln_gamma_s_MR=numpy.zeros((1,org_molecules),)
    #for molecule_step in range(org_molecules):
    #    for k in range(max_group_num_org_main):
    #        Ln_gamma_s_MR[0,molecule_step]=Ln_gamma_s_MR[0,molecule_step]+org_group_stoich_new_main[molecule_step,k]*ln_gamma_k_MR[0,k]

    #Now sum all of the individual functional groups together to get the final activity coefficient
    Ln_gamma_s_MR = {
        compound: sum(ln_gamma_k_MR[group]*m_org[compound][group] for group in non_zero_groups)
    for compound, abundance in organic_compounds.items()
    }

    #Activity coefficient for ions
    #Generic calculations
    summation1_ion = {
    ion: sum(Bki[group][ion]*mole_frac_non_zero_groups[group]*(1.0/M_solv_mix)
    for group in non_zero_groups)
        for ion, value in ion_molalities.items()
    }

    pre_summation2_ion = sum(Bki_dash[group][ion]*mole_frac_non_zero_groups[group]*value
    for ion, value in ion_molalities.items()
    for group in non_zero_groups)

    summation2_ion = {
    ion:pre_summation2_ion*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)/2.0*M_solv_mix)
    for ion, value in ion_molalities.items()
    }

    pre_summation4_ion = sum(
    Bca_dash[cation][anion]*c_val*a_val
    for anion, a_val in anion_molality.items()
    for cation, c_val in cation_molality.items()
    )
    summation4_ion = {
    ion: pre_summation4_ion*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)*0.5)
    for ion, value in ion_molalities.items()
    }
    pre_summation5_ion=sum(abs(data.AIOMFAC_ION_CHARGE[ion])*value for ion, value in ion_molalities.items())
    pre_summation6_ion = pre_summation5_ion

    summation6_ion = {
        ion:sum((Cca[cation][anion]*abs(data.AIOMFAC_ION_CHARGE[ion])+Cca_dash[cation][anion]*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)*0.5)*pre_summation6_ion)*c_val*a_val
    for anion, a_val in anion_molality.items()
    for cation, c_val in cation_molality.items())
    for ion, value in ion_molalities.items()
    }

    #Cation/anion specific calculations
    summation3_ion_cation = {
    cation: sum (Bca[cation][anion]*value for anion, value in anion_molality.items())
    for cation, value_c in cation_molality.items()
    }
    summation3_ion_anion = {
    anion: sum (Bca[cation][anion]*value for cation, value in cation_molality.items())
    for anion, value_a in anion_molality.items()
    }
    summation3_ion=summation3_ion_cation.copy()
    summation3_ion.update(summation3_ion_anion)

    summation5_ion_cation = {
        cation : sum(Cca[cation][anion]*value*pre_summation5_ion for anion, value in anion_molality.items())
    for cation, value_c in cation_molality.items()
    }
    summation5_ion_anion = {
    anion : sum(Cca[cation][anion]*value*pre_summation5_ion for cation, value in cation_molality.items())
    for anion, value_a in anion_molality.items()
    }
    summation5_ion=summation5_ion_cation.copy()
    summation5_ion.update(summation5_ion_anion)

    summation7_ion_cation_1 = {
        cation: sum(Rc_cdash[cation][ion]*value for ion, value in cation_molality.items())
    for cation, value_c in cation_molality.items()
    }
    summation7_ion_cation_2 = {
    cation: sum(Qc_cdash[cation][c][a]*c_val*a_val
    for c,c_val in cation_molality.items()
    for a,a_val in anion_molality.items())
    for cation, value_c in cation_molality.items()
    }
    summation7_ion_cation = {
    cation: summation7_ion_cation_1[cation]+summation7_ion_cation_2[cation]
    for cation, value_c in cation_molality.items()
    }
    summation7_ion_anion = {
    anion: sum(Qc_cdash[c][cdash][anion]*c_val*cdash_val
    for c, c_val in cation_molality.items()
    for cdash, cdash_val in cation_molality.items())
    for anion , a_val in anion_molality.items()
    }
    summation7_ion=summation7_ion_cation.copy()
    summation7_ion.update(summation7_ion_anion)

    Ln_gamma_i_MR={
    ion: summation1_ion[ion]+summation2_ion[ion]+summation3_ion[ion]+summation4_ion[ion]+\
        summation5_ion[ion]+summation6_ion[ion]+summation7_ion[ion]
    for ion, value in ion_molalities.items()
        if value > 0
    }

    #--------------------------------------------------------------------------------------------
    #Now calculate the LR portion of the activity coefficient for both organics and ions
    eo=1.602177e-19; #elementary charge (C)
    N_A=6.02213e23; #Avogadro's number
    k=1.381e-23; #Boltzmann constant (J/K)
    a=10e-10; #closest approach parameter (m)
    epsilon_o=8.854187817e-12; #permittivitty of a vacuum (C^2/Jm)
    epsilon_wr=81; #relative permittivity water
    D=78.54 #Dielecric constant for water - taken directly from the AIOMFAC code
    #Note this changes with temperature
    dens_w=997; #density of water (kg/m3) at 298K
    #A=1.327757e5*((numpy.power(dens_w,0.5))/(numpy.power((epsilon_wr*T),1.5)));
    A=1.327757e5*((pow(dens_w,0.5))/(pow((D*temperature),1.5)));
    #print'A (LR)',A
    #b=6.359696*(numpy.power(dens_w,0.5))/(numpy.power((epsilon_wr*T),0.5));
    b=6.359696*(pow(dens_w,0.5))/(pow((D*temperature),0.5));

    summation_org=(1+b*pow(Ionic_strength,0.5)-(1/(1+b*pow(Ionic_strength,0.5)))-2.0*log(1+b*pow(Ionic_strength,0.5)));
    summation_ion_1=A*pow(Ionic_strength,0.5);
    summation_ion_2=1+b*pow(Ionic_strength,0.5);


    #--Organic solvent activity coefficients---
    Ln_gamma_s_LR={
        compound: ((2.0*A*(compound.molwt*1e-3))/(pow(b,3.0)))*summation_org
        for compound, abundance in organic_compounds.items()
        }

    #--Ionic activity coefficients--
    Ln_gamma_i_LR={
        ion: (-1.0*pow((abs(data.AIOMFAC_ION_CHARGE[ion])),2.0)*summation_ion_1)/summation_ion_2 #on the mole fraction scale
        for ion, value in ion_molalities.items()
        }

    #Now calculate the conversion factor to change the reference state from mole fraction to molality
    total_abundance_org = sum(value for key, value in organic_compounds.iteritems())
    x_i_org = {
        compound: value/total_abundance_org for compound, value in organic_compounds.items()
        }
    mean_mw_solvent=sum(x_i_org[compound]*(compound.molwt*1e-3) for compound, value in organic_compounds.items())
    sum_molalities=sum(value for ion, value in ion_molalities.items())
    ionic_conversion_factor=log((0.01801528/mean_mw_solvent)+0.01801528*sum_molalities)

    Ln_gamma_i_MR_LR= {
        ion: Ln_gamma_i_MR[ion]+Ln_gamma_i_LR[ion]-ionic_conversion_factor
        for ion, value in ion_molalities.items()
        if value > 0
        }
    #Now map this ontp pybel objects as keys rather than integers

    Ln_gamma_i_MR_LR_keys = {
        ion: Ln_gamma_i_MR_LR[old_key] #q_k_i[compound].get(group1, 0)
        for ion, matches in m_inorg.items()
        for old_key, count in matches.items()
        if count >0
        }
    Ln_gamma_s_MR_LR = {
        compound: Ln_gamma_s_LR[compound]+Ln_gamma_s_MR[compound]
        for compound, abundance in organic_compounds.items()
        }
        
    # Create a dictionary that holds both ionic and organic values
    
    Ln_gamma_tot_MR_LR={}
    for compound, abundance in organic_compounds.items():
        Ln_gamma_tot_MR_LR[compound] = Ln_gamma_s_MR_LR[compound]
    
    for ion, matches in m_inorg.items():
        Ln_gamma_tot_MR_LR[ion]= Ln_gamma_i_MR_LR_keys[ion]


    #NOTE we only return the organic activity coefficients here. For including inorganic
    #ions, this will also need to be passed
    return Ln_gamma_tot_MR_LR




def calculate_activities_sr(organic_compounds, inorganic_ions, temperature):
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
         #molar_mass=compound.molwt
         )
      for compound, abundance in organic_compounds.items()
      ]
    m.extend([generic_class(
         compound=compound,
         abundance=abundance
         #molar_mass=compound.molwt
         )
      for compound, abundance in inorganic_ions.items()
      ])

    #BCoa = sum(c.abundance for c in m)
    salts = aiomfac_salts(inorganic_ions)
    core_ion = sum(inorganic_ions.values())

    #Now calculate activity coefficients associated with each component
    #return it as a value associated with each

    Activity_coefficients_sr=aiomfac_sr(organic_compounds, inorganic_ions, temperature)
    step=0
    for c in m:
        x=m[step].compound
        c.update(Activity_coefficients_sr[x])
        step+=1

    #for compound, abundance in m.items(0

    return Activity_coefficients_sr, m

def calculate_activities_full(organic_compounds, inorganic_ions, temperature):
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
         #molar_mass=compound.molwt
         )
      for compound, abundance in organic_compounds.items()
      ]
    m.extend([generic_class(
         compound=compound,
         abundance=abundance
         #molar_mass=compound.molwt
         )
      for compound, abundance in inorganic_ions.items()
      ])

    #BCoa = sum(c.abundance for c in m)
    salts = aiomfac_salts(inorganic_ions)
    core_ion = sum(inorganic_ions.values())

    #Now calculate activity coefficients associated with each component
    #Note here we combine short range (SR) with long range+mid range (MR_LR)
    #If you want to call the activity module seperately, you need to bare this in
    #mind and always combine the two. You can see how this is done in partition_models.partition_model

    Activity_coefficients_sr=aiomfac_sr(organic_compounds, inorganic_ions, temperature)
    Activity_coefficients_mr_lr=aiomfac_mr(organic_compounds, inorganic_ions, temperature)
    #pdb.set_trace()
    step=0
    for c in m:
        x=m[step].compound
        c.update(exp(Activity_coefficients_sr[x]+Activity_coefficients_mr_lr[x]))
        step+=1

    #for compound, abundance in m.items(0

    return Activity_coefficients_mr_lr, m

def calculate_activities_org(organic_compounds, temperature):
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
         #molar_mass=compound.molwt
         )
      for compound, abundance in organic_compounds.items()
      ]
    #m.extend([generic_class(
    #     compound=compound,
    #     abundance=abundance
         #molar_mass=compound.molwt
    #     )
    #  for compound, abundance in inorganic_ions.items()
    #  ])

    #BCoa = sum(c.abundance for c in m)
    #salts = aiomfac_salts(inorganic_ions)
    #core_ion = sum(inorganic_ions.values())

    #Now calculate activity coefficients associated with each component
    #return it as a value associated with each

    Activity_coefficients_sr=aiomfac_sr(organic_compounds, {}, temperature)
    #Activity_coefficients_mr_lr=aiomfac_mr(organic_compounds, inorganic_ions, temperature)
    step=0
    for c in m:
        x=m[step].compound
        c.update(exp(Activity_coefficients_sr[x]))
        step+=1

    #for compound, abundance in m.items(0

    return Activity_coefficients_sr, m



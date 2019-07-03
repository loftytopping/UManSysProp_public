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
from . import groups as groups_func
from .forms import smiles
from functools import wraps

#-------------------------------------------------
#Development only - not kept in standard release
import pdb
#import cPickle
import os.path
import uuid
import weakref
import time
import numpy as np
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
        compound: groups_func.aiomfac(compound)
        for compound in compounds
        }


def aiomfac_inorganic(compounds):
    return {
        compound: groups_func.matches(data.AIOMFAC_ION_SMARTS, compound)
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

#def memo(func):
#    cache = {}
#    @wraps(func)
#    def wrap(*args):
#        if args not in cache:
#            pdb.set_trace()
#            cache[args] = func(*args)
#        return cache[args]
#    return wrap


#The following functions are designed to work with  'fast' version of aiomfac
#that might be used in box models 


def aiomfac_sr_persistent(organic_compounds, inorganic_ions, temperature,species_dict2array,Pybel_object_to_name):
    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    m_org = aiomfac_organic(organic_compounds)
    m_inorg = aiomfac_inorganic(inorganic_ions)
    
    # There can be a disconnect between the number of compounds in 'm' and those in species_dict2array
    # The former deals with ensuring the matrix size of components ism correct. Thus, this 
    # dictates the size of created arrays
    
    number_compounds = max(species_dict2array.values())+1

    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    
    non_zero_groups = {
            group
            for compound, matches in m.items()
            for group, count in matches.items()
            if count > 0
            }
    #This array tells us which column corresponds to which group
    non_zero_groups_list = []
    for group in non_zero_groups:
        if group not in non_zero_groups_list:
            non_zero_groups_list.append(group)
    
    q_k_i = {
            compound: {
                group: count * data.AIOMFAC_QI[group]
                for group, count in matches.items()
                }
            for compound, matches in m.items()
            }  
    q_k_i_np=np.zeros((number_compounds,len(non_zero_groups)),)
    #pdb.set_trace()
    for compound, matches in m.items():
        for group, count in matches.items():
            if count > 0:
                try:
                    q_k_i_np[species_dict2array[Pybel_object_to_name[compound]],non_zero_groups_list.index(group)]=q_k_i[compound][group]
                except:
                    pdb.set_trace()
            
    q_i = {
        compound: groups_func.aggregate_matches(matches, data.AIOMFAC_QI)
        for compound, matches in m.items()
        }
    q_i_np = np.zeros(number_compounds,)
    for compound in q_i.keys():
        q_i_np[species_dict2array[Pybel_object_to_name[compound]]]=q_i[compound]
        

    r_i = {
        compound: groups_func.aggregate_matches(matches, data.AIOMFAC_RI)
        for compound, matches in m.items()
        }
    r_i_np = np.zeros(number_compounds,)
    for compound in r_i.keys():
        r_i_np[species_dict2array[Pybel_object_to_name[compound]]]=r_i[compound]

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

    #This is also passed out
    #pdb.set_trace()
    u_i_m_n_np=np.zeros((number_compounds,len(non_zero_groups),len(non_zero_groups)),)
    for compound, matches in m.items():
        for group1 in non_zero_groups:
            for group2 in non_zero_groups:
                u_i_m_n_np[species_dict2array[Pybel_object_to_name[compound]],non_zero_groups_list.index(group1),non_zero_groups_list.index(group2)]=u_i_m_n[compound][group1][group2]     
    #pdb.set_trace()

    s_i_n = {
        compound: {
            group1: sum(t_i_m_n[compound][group2][group1] for group2 in non_zero_groups)
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }
    s_i_n_np=np.zeros((number_compounds,len(non_zero_groups)),)
    for compound, matches in m.items():
        for group1 in non_zero_groups:
            s_i_n_np[species_dict2array[Pybel_object_to_name[compound]],non_zero_groups_list.index(group1)]=s_i_n[compound][group1]
    #pdb.set_trace()
    brac1 = {
        compound: {
            group1: sum(u_i_m_n[compound][group2][group1]/s_i_n[compound][group2] for group2 in non_zero_groups)
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }
    brac1_np=np.zeros((number_compounds,len(non_zero_groups)),)
    for compound, matches in m.items():
        for group1 in non_zero_groups:
            brac1_np[species_dict2array[Pybel_object_to_name[compound]],non_zero_groups_list.index(group1)]=brac1[compound][group1]
    #pdb.set_trace()    
    
    persistent_dict=dict()
    persistent_dict['non_zero_groups_list']=non_zero_groups_list
    persistent_dict['num_non_zero_groups']=len(non_zero_groups)
    persistent_dict['q_k_i_np']=q_k_i_np
    persistent_dict['q_i_np']=q_i_np
    persistent_dict['r_i_np']=r_i_np
    persistent_dict['u_i_m_n_np']=u_i_m_n_np
    persistent_dict['s_i_n_np']=s_i_n_np
    persistent_dict['brac1_np']=brac1_np
    
    return persistent_dict
            
def aiomfac_sr_quick1(abundance,persistent_dict, temperature):    
    
    non_zero_groups_list=persistent_dict['non_zero_groups_list']
    q_k_i_np=persistent_dict['q_k_i_np']
    q_i_np=persistent_dict['q_i_np']
    r_i_np=persistent_dict['r_i_np']
    u_i_m_n_np=persistent_dict['u_i_m_n_np']
    s_i_n_np=persistent_dict['s_i_n_np']
    brac1=persistent_dict['brac1_np']
    
    numb_nonzero_groups=len(non_zero_groups_list)
    numb_compounds=len(abundance)

    x_i=abundance/(np.sum(abundance))

    #pdb.set_trace()

    #Q_old = sum(x_i[compound]*q_i[compound] for compound, matches in m.items())
    Q=np.sum(x_i.flatten()*q_i_np.flatten())

    #omega_i_old = {
    #    compound: log(q_i[compound]/Q) for compound, matches  in m.items()
    #    }
    omega_i=np.log(q_i_np/Q)

    #R_old = sum(x_i[compound]*r_i[compound] for compound, matches  in m.items())
    R=np.sum(x_i.flatten()*r_i_np.flatten())

    #row_i_old = {
    #    compound: r_i[compound]/R for compound, matches  in m.items()
    #    }
    row_i=r_i_np/R

    #P_i_old = {
    #    compound: log(r_i[compound]/R) for compound, matches  in m.items()
    #    }
    P_i = np.log(r_i_np/R)

    #delta_i_old = {
    #    compound: row_i[compound]*(5.0*Q-1.0) for compound, matches  in m.items()
    #    }
    delta_i = np.multiply(row_i,(5.0*Q-1.0))

    #cross_i_old = {
    #    compound: 5.0*q_i[compound] for compound, matches  in m.items()
    #    }
    cross_i = np.multiply(q_i_np,5.0)

    #Ln_gamma_i_C_old = {
    #    compound: 1.0+delta_i[compound]+P_i[compound]+cross_i[compound]*(omega_i[compound]-P_i[compound]-1.0) for compound, matches  in m.items()
    #    }
    Ln_gamma_i_C = 1.0+delta_i+P_i+np.multiply(cross_i,(omega_i-P_i-1.0))
    
    #pdb.set_trace()
    #Now do the residual calculation [ignoring ions as seperate outputs for now!]
    #uu_m_n = {
    #   group1: {
    #       group2: sum(x_i[compound]*u_i_m_n[compound.title][group1][group2] for compound, matches  in m.items())
    #       for group2 in non_zero_groups
    #       }
    #    for group1 in non_zero_groups
    #    }
    #pdb.set_trace()
    uu_m_n=np.zeros((numb_nonzero_groups,numb_nonzero_groups),)
    for grp1 in range(numb_nonzero_groups):
        for grp2 in range(numb_nonzero_groups):
            uu_m_n[grp1,grp2]=np.sum(x_i[:,0]*u_i_m_n_np[:,grp1,grp2])

    #ss_n = {
    #    group1: sum(x_i[compound]*s_i_n[compound][group1] for compound, matches in m.items())
    #    for group1 in non_zero_groups
    #    }
    #ss_n=np.zeros((numb_nonzero_groups,1),)
    ss_n=np.sum(x_i*s_i_n_np,axis=0)

    #xx_weird_i_n = {
    #    compound: {
    #        group1: log(s_i_n[compound][group1]/ss_n[group1])
    #        for group1 in non_zero_groups
    #        }
    #    for compound, matches in m.items()
    #    }
    #xx_weird_i_n=np.zeros((numb_compounds,numb_nonzero_groups),)
    
    xx_weird_i_n=np.log(s_i_n_np/ss_n)

    #brac2 = {
    #    compound: {
    #        group1: sum(uu_m_n[group2][group1]/ss_n[group2] for group2 in non_zero_groups)
    #        for group1 in non_zero_groups
    #        }
    #    for compound, matches in m.items()
    #    }
    #brac2 = np.zeros((numb_compounds,numb_nonzero_groups),)
    temp=np.sum((uu_m_n.transpose()/ss_n),axis=1).transpose()
    brac2 = np.tile(temp,(numb_compounds,1))

    #pdb.set_trace()
    
    #summation2 = {
    #    compound.title: {
    #        group1: xx_weird_i_n[compound][group1]-omega_i[compound]+brac1[compound.title][group1]-brac2[compound][group1]
    #         for group1 in non_zero_groups
    #        }
    #    for compound, matches in m.items()
    #    }
    # We dont need to loop through all items.
    summation2 = xx_weird_i_n-omega_i[:,None]+brac1-brac2
        
    #compared to the original code we have taken out the use of 'summation1'
    #Ln_gamma_i_R = {
    #    compound: sum(q_k_i[compound].get(group1, 0)*summation2[compound][group1] for group1 in# non_zero_groups)
    #    for compound, matches in m.items()
    #    }
    #pdb.set_trace()
    Ln_gamma_i_R = np.sum(np.prod([q_k_i_np,summation2],axis=0),axis=1)

    Ln_gamma_i_SR = Ln_gamma_i_R + Ln_gamma_i_C

    return Ln_gamma_i_SR

def aiomfac_sr_quick(abundance,persistent_dict, temperature):    

    non_zero_groups_list=persistent_dict['non_zero_groups_list']
    q_k_i_np=persistent_dict['q_k_i_np']
    q_i_np=persistent_dict['q_i_np']
    r_i_np=persistent_dict['r_i_np']
    u_i_m_n_np=persistent_dict['u_i_m_n_np']
    s_i_n_np=persistent_dict['s_i_n_np']
    brac1=persistent_dict['brac1_np']
    num_orgs=persistent_dict['num_orgs']
    num_ions=persistent_dict['num_ions']
    molw_array_org=persistent_dict['molw_array_org']
    
    numb_nonzero_groups=len(non_zero_groups_list)
    numb_compounds=len(abundance)    
    
    #total_abundance = sum(value for key, value in list_compounds.iteritems())

    x_i=abundance/(np.sum(abundance))

    Q=np.sum(x_i.flatten()*q_i_np.flatten())
     
    omega_i=np.log(q_i_np/Q)

    R=np.sum(x_i.flatten()*r_i_np.flatten())

    row_i=r_i_np/R

    P_i = np.log(r_i_np/R)

    delta_i = np.multiply(row_i,(5.0*Q-1.0))

    cross_i = np.multiply(q_i_np,5.0)

    Ln_gamma_i_C = 1.0+delta_i+P_i+np.multiply(cross_i,(omega_i-P_i-1.0))
    
    uu_m_n=np.zeros((numb_nonzero_groups,numb_nonzero_groups),)
    for grp1 in range(numb_nonzero_groups):
        for grp2 in range(numb_nonzero_groups):
            uu_m_n[grp1,grp2]=np.sum(x_i[:,0]*u_i_m_n_np[:,grp1,grp2])

    ss_n=np.sum(x_i*s_i_n_np,axis=0)

    xx_weird_i_n=np.log(s_i_n_np/ss_n)

    brac2 = np.zeros((numb_compounds,numb_nonzero_groups),)
    temp=np.sum((uu_m_n.transpose()/ss_n),axis=1).transpose()
    brac2 = np.tile(temp,(numb_compounds,1))

    summation2 = xx_weird_i_n-omega_i[:,None]+brac1-brac2

    Ln_gamma_i_R = np.sum(np.prod([q_k_i_np,summation2],axis=0),axis=1)

    Ln_gamma_i_SR = Ln_gamma_i_R + Ln_gamma_i_C
    
    if num_ions > 0:
        ion_abundance=abundance[num_orgs::]
        solvent_kg=np.sum(molw_array_org*np.array(abundance[0:num_orgs,0]))*1.0e-3
        ion_molalities = ion_abundance / solvent_kg
        x_i_org_persistent=np.array(abundance[0:num_orgs,0])/(np.sum(np.array(abundance[0:num_orgs,0])))
        #pdb.set_trace()
        mean_mw_solvent_persistent=np.sum(np.multiply(x_i_org_persistent,molw_array_org)*1e-3)  
        #pdb.set_trace()
        sum_molalities_persistent=np.sum(ion_molalities) 
        #pdb.set_trace()
        ionic_conversion_factor=log((0.01801528/mean_mw_solvent_persistent)+0.01801528*sum_molalities_persistent)
 
        Ln_gamma_i_SR[num_orgs:num_orgs+num_ions]=Ln_gamma_i_SR[num_orgs:num_orgs+num_ions]-ionic_conversion_factor

    return Ln_gamma_i_SR
    


def aiomfac_sr_quick_test(abundance,persistent_dict, temperature,organic_compounds, inorganic_ions):    
    
    non_zero_groups_list=persistent_dict['non_zero_groups_list']
    q_k_i_np=persistent_dict['q_k_i_np']
    q_i_np=persistent_dict['q_i_np']
    r_i_np=persistent_dict['r_i_np']
    u_i_m_n_np=persistent_dict['u_i_m_n_np']
    s_i_n_np=persistent_dict['s_i_n_np']
    brac1=persistent_dict['brac1_np']
    
    numb_nonzero_groups=len(non_zero_groups_list)
    numb_compounds=len(abundance)

    #m = aiomfac_organic(organic_compounds)
    #m.update(aiomfac_inorganic(inorganic_ions))

    #m_org = aiomfac_organic(organic_compounds)
    #m_inorg = aiomfac_inorganic(inorganic_ions)

    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    # XXX What do any of these variables mean?! Unfortunately the former coder
    # left no clue beyond their name ...

    # XXX Need to exclude any groups which have zero matches in the following
    # calculation; the original does this as an "optimization" and yet it
    # produces a different answer. Hmmm ...
    
    #Check to see if pickled unifac objects already exist
    
    #pdb.set_trace()
    
    
    #non_zero_groups_old = {
    #    group
    #    for compound, matches in m.items()
    #    for group, count in matches.items()
    #    if count > 0
    #    }

    #q_k_i_old = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count in matches.items()
    #        }
    #    for compound, matches in m.items()
    #    }  
        
    #q_k_i = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count in matches.items()
    #        #for group, count in non_zero_groups
    #        }
    #    for compound, matches in m.items()
    #    }

    #q_k_i = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count

    #q_i_old = {
    #    compound: groups.aggregate_matches(matches, data.AIOMFAC_QI)
    #    for compound, matches in m.items()
    #    }

    #r_i_old = {
    #    compound: groups.aggregate_matches(matches, data.AIOMFAC_RI)
    #    for compound, matches in m.items()
    #    }

    #t_i_m_n_old = {
    #    compound: {
    #        group1: {
    ##            group2: q_k_i_old[compound].get(group1, 0) * exp(
    #                -data.AIOMFAC_SR_INTERACTIONS[main_group1][main_group2] /
    #                temperature)
    #            for group2 in non_zero_groups_old
    #            for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
    #            }
    #        for group1 in non_zero_groups_old
    #        for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
    #        }
    #    for compound, matches in m.items()
    #    }

    #u_i_m_n_old = {
    #    compound: {
    #        group1: {
    ##            group2: q_k_i_old[compound].get(group1, 0) * exp(
    #                -data.AIOMFAC_SR_INTERACTIONS[main_group2][main_group1] /
    #                temperature)
    #            for group2 in non_zero_groups_old
    #            for main_group2 in (data.AIOMFAC_MAIN_GROUP[group2],)
    #            }
    #        for group1 in non_zero_groups_old
    #        for main_group1 in (data.AIOMFAC_MAIN_GROUP[group1],)
    #        }
    #    for compound, matches in m.items()
    #    }

    #s_i_n_old = {
    #    compound: {
    #        group1: sum(t_i_m_n_old[compound][group2][group1] for group2 in non_zero_groups_old)
    #        for group1 in non_zero_groups_old
    #        }
    #    for compound, matches in m.items()
    #    }
    
    #brac1_old = {
    #    compound: {
    #        group1: sum(u_i_m_n_old[compound][group2][group1]/s_i_n_old[compound][group2] for group2 in non_zero_groups_old)
    #        for group1 in non_zero_groups_old
    #        }
    #    for compound, matches in m.items()
    #    }
                       
            
    #pdb.set_trace()
            
    #water = [c for c in organic_compounds if str(c).strip() == 'O'][0]

    #Now adding the rest of the calculation to complete predictions of activity
    #coefficients from the UNIFAC portion of AIOMFAC [sr]

    ##mole fraction
    #list_compounds=organic_compounds.update(inorganic_ions)
    
    
    total_abundance = sum(value for key, value in list_compounds.iteritems())
    #x_i_old = {
    #    compound: value/total_abundance for compound, value in list_compounds.iteritems()
    #    }


    x_i=abundance/(np.sum(abundance))

    #pdb.set_trace()

    #Q_old = sum(x_i_old[compound]*q_i_old[compound] for compound, matches in m.items())
    Q=np.sum(x_i.flatten()*q_i_np.flatten())
     
    #omega_i_old = {
    #    compound: log(q_i_old[compound]/Q_old) for compound, matches  in m.items()
    #    }
    omega_i=np.log(q_i_np/Q)

    #R_old = sum(x_i_old[compound]*r_i_old[compound] for compound, matches  in m.items())
    R=np.sum(x_i.flatten()*r_i_np.flatten())

    #row_i_old = {
    #    compound: r_i_old[compound]/R_old for compound, matches  in m.items()
    #    }
    row_i=r_i_np/R

    #P_i_old = {
    #    compound: log(r_i_old[compound]/R_old) for compound, matches  in m.items()
    #    }
    P_i = np.log(r_i_np/R)

    #delta_i_old = {
    #    compound: row_i_old[compound]*(5.0*Q_old-1.0) for compound, matches  in m.items()
    #    }
    delta_i = np.multiply(row_i,(5.0*Q-1.0))

    #cross_i_old = {
    #    compound: 5.0*q_i_old[compound] for compound, matches  in m.items()
    #    }
    cross_i = np.multiply(q_i_np,5.0)

    #Ln_gamma_i_C_old = {
    #   compound: 1.0+delta_i_old[compound]+P_i_old[compound]+cross_i_old[compound]*(omega_i_old[compound]-P_i_old[compound]-1.0) for compound, matches # in m.items()
    #        }
    Ln_gamma_i_C = 1.0+delta_i+P_i+np.multiply(cross_i,(omega_i-P_i-1.0))
    
    #pdb.set_trace()
    #Now do the residual calculation [ignoring ions as seperate outputs for now!]
    #uu_m_n_old = {
    #   group1: {
    #       group2: sum(x_i_old[compound]*u_i_m_n_old[compound][group1][group2] for compound, matches  in m.items())
    #       for group2 in non_zero_groups_old
    #       }
    #    for group1 in non_zero_groups_old
    #    }
    #pdb.set_trace()
    uu_m_n=np.zeros((numb_nonzero_groups,numb_nonzero_groups),)
    for grp1 in range(numb_nonzero_groups):
        for grp2 in range(numb_nonzero_groups):
            uu_m_n[grp1,grp2]=np.sum(x_i[:,0]*u_i_m_n_np[:,grp1,grp2])

    #pdb.set_trace()
    #ss_n_old = {
    #    group1: sum(x_i_old[compound]*s_i_n_old[compound][group1] for compound, matches in m.items())
    #    for group1 in non_zero_groups_old
    #    }
    #ss_n=np.zeros((numb_nonzero_groups,1),)
    ss_n=np.sum(x_i*s_i_n_np,axis=0)

    

    #xx_weird_i_n_old = {
    #    compound: {
    #        group1: log(s_i_n_old[compound][group1]/ss_n_old[group1])
    #        for group1 in non_zero_groups_old
    #        }
    #    for compound, matches in m.items()
    #    }
    #xx_weird_i_n=np.zeros((numb_compounds,numb_nonzero_groups),)
    
    xx_weird_i_n=np.log(s_i_n_np/ss_n)

    #pdb.set_trace()

    #brac2_old = {
    #    compound: {
    #        group1: sum(uu_m_n_old[group2][group1]/ss_n_old[group2] for group2 in non_zero_groups_old)
    #        for group1 in non_zero_groups_old
    #        }
    #    for compound, matches in m.items()
    #    }
    brac2 = np.zeros((numb_compounds,numb_nonzero_groups),)
    temp=np.sum((uu_m_n.transpose()/ss_n),axis=1).transpose()
    brac2 = np.tile(temp,(numb_compounds,1))

    #pdb.set_trace()
    
    #summation2_old = {
    #    compound: {
    #        group1: xx_weird_i_n_old[compound][group1]-omega_i_old[compound]+brac1_old[compound][group1]-brac2_old[compound][group1]
    #         for group1 in non_zero_groups_old
    #        }
    #    for compound, matches in m.items()
    #    }
    # We dont need to loop through all items.
    summation2 = xx_weird_i_n-omega_i[:,None]+brac1-brac2

    #pdb.set_trace()
        
    #compared to the original code we have taken out the use of 'summation1'
    #Ln_gamma_i_R_old = {
    #    compound: sum(q_k_i_old[compound].get(group1, 0)*summation2_old[compound][group1] for group1 in non_zero_groups_old)
    #    for compound, matches in m.items()
    #    }
    #pdb.set_trace()
    Ln_gamma_i_R = np.sum(np.prod([q_k_i_np,summation2],axis=0),axis=1)

    #pdb.set_trace()

    #final activity coefficient [without correction for ions]
    Ln_gamma_i_SR = Ln_gamma_i_R + Ln_gamma_i_C

    return Ln_gamma_i_SR

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
    
    #Check to see if pickled unifac objects already exist
    
    #pdb.set_trace()
    
    
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
        
    #q_k_i = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count in matches.items()
    #        #for group, count in non_zero_groups
    #        }
    #    for compound, matches in m.items()
    #    }

    #q_k_i = {
    #    compound: {
    #        group: count * data.AIOMFAC_QI[group]
    #        for group, count

    q_i = {
        compound: groups_func.aggregate_matches(matches, data.AIOMFAC_QI)
        for compound, matches in m.items()
        }

    r_i = {
        compound: groups_func.aggregate_matches(matches, data.AIOMFAC_RI)
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
        for compound, matches in m.items()
        }
    
    brac1 = {
        compound: {
            group1: sum(u_i_m_n[compound][group2][group1]/s_i_n[compound][group2] for group2 in non_zero_groups)
            for group1 in non_zero_groups
            }
        for compound, matches in m.items()
        }
                       
            
        #pdb.set_trace()
            
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
    m_ions = groups_func.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
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
                compound.title: value/total_abundance_org for compound, value in organic_compounds.items()
                }
            mean_mw_solvent=sum(x_i_org[compound.title]*(compound.molwt*1e-3) for compound, value in organic_compounds.items())
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
                
    #pdb.set_trace()
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

    #total_non_zero_groups=0
    #conc_non_zero_groups={}
    #for compound, matches in m_org.items():
    #    for group, count in matches.items():
    #        if count > 0:
    #            #conc_non_zero_groups[group]=conc_non_zero_groups[group]+count
    #            total_non_zero_groups=total_non_zero_groups+count

    #only over the organic functional groups
    conc_non_zero_groups = {
        groups: sum(m_org[compound].get(groups, 0)*organic_compounds[compound] for compound, matches in m_org.items())
        for groups in non_zero_groups
                           }
                           
    total_non_zero_groups=sum(conc_non_zero_groups.values())                       
    #again, only over the organic compounds
    mole_frac_non_zero_groups = {
        groups: conc_non_zero_groups[groups]/total_non_zero_groups
        for groups in non_zero_groups
                                }
    #now calculate properties associated with the ions - molalities
    #pdb.set_trace()
    m_ions = groups_func.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
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
    #pdb.set_trace()
    return Ln_gamma_tot_MR_LR
    
def aiomfac_mr_persistent(organic_compounds, inorganic_ions, temperature,species_dict2array,Pybel_object_to_name,ion_dict2array,num_species,abundance_array,cation_index,anion_index):
    
    #The purpose of this function is the create the persistant variables that are used throughout the calculation of MR terms
    #It uses the flexible dictionary nature of the original development but provides variables for use with the 'fast version'

    #Also adding the option to write the contents of a Fortran equivalent sub-module for use
    #in the box model. Creat a module that just defines values if integers used to create
    #arrays and matrices in proceeding modules created by the box-model    
    persistant_data=dict()
    
    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    # There can be a disconnect between the number of compounds in 'm' and those in species_dict2array
    # The former deals with ensuring the matrix size of components ism correct. Thus, this 
    # dictates the size of created arrays
    
    number_compounds = max(species_dict2array.values())+1

    m_org = aiomfac_organic(organic_compounds)
    m_inorg = aiomfac_inorganic(inorganic_ions)
    #now calculate properties associated with the ions - molalities
    m_ions = groups_func.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
    
    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    #Mole fraction of all compounds

    #For the persistent version we need some way of extracting ions from the abundance array
    #We will have to do this from the calling script.
    #pdb.set_trace()
    ion_abundance=abundance_array[num_species::]

    non_zero_groups_old = {
        group
        for compound, matches in m_org.items()
        for group, count in matches.items()
        if count > 0
        }
    #--Persistent variable--
    non_zero_groups_list = []
    for group in non_zero_groups_old:
        if group not in non_zero_groups_list:
            non_zero_groups_list.append(group)

    persistant_data['non_zero_groups_list']=non_zero_groups_list

    #pdb.set_trace()

    #Now create a matrix that simply records the flag on a nonzero group or not
    #Persistent variable
    non_zero_groups_flag_old = {
        compound : { group : count
        for group, count in matches.items()
        }
        for compound, matches in m_org.items()
    }
    #--Persistent variable--
    non_zero_groups_flag = np.zeros((num_species,len(non_zero_groups_list)),)
    #pdb.set_trace()
    for compound, matches in m_org.items():
        #pdb.set_trace()
        for group, count in matches.items():
            #pdb.set_trace()
            if group in non_zero_groups_list:
                non_zero_groups_flag[species_dict2array[Pybel_object_to_name[compound]],non_zero_groups_list.index(group)]=count
                 
    persistant_data['non_zero_groups_flag']=non_zero_groups_flag

    #--Persistent variable--
    #Persistent variable - [molw_array]
    num_orgs=num_species
    #--Persistent variable--
    molw_array_org=np.zeros((num_orgs),)
    #pdb.set_trace()
    for key, value in organic_compounds.items():
        #pdb.set_trace()
        molw_array_org[species_dict2array[Pybel_object_to_name[key]]]=key.molwt
        #print(Pybel_object_to_name[key])
        #print(key.molwt)
        #if species_dict2array[Pybel_object_to_name[key]] == 305:
    #pdb.set_trace()
        #if key.molwt == 0.0:
        #   pdb.set_trace()
        
    persistant_data['num_orgs']=num_orgs
    persistant_data['num_non_zero_groups']=len(non_zero_groups_list)
    persistant_data['molw_array_org']=molw_array_org    
    
    #Cation/Anion molality
    cations = {
        ion
        for ion, count in m_ions.items() #ion is a numeric key, not a pybel object
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    cations=list(cations)
    #Now setup a dictionary that creates cation array indices
    step=0
    cation2array = dict()
    for ion in cations:
        cation2array[ion]=step
        step+=1
    
    persistant_data['cation2array']=cation2array   

    anions  = {
        ion
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    anions=list(anions)
    step=0
    anion2array = dict()
    for ion in anions:
        anion2array[ion]=step
        step+=1
    
    persistant_data['anion2array']=anion2array  

    #Persistent variable
    ion_charge=np.zeros((len(ion_abundance),1),)
    for ion, info in m_inorg.items():
        for key, value in info.items():
            if value > 0:
                key_rec=key
                ion_charge[ion_dict2array[ion],0]=data.AIOMFAC_ION_CHARGE[key_rec]
    
    persistant_data['ion_charge']=ion_charge

    #pdb.set_trace()    
    #Persistent version
    
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
        for groups in non_zero_groups_old
        }

    b2_ki_temp = {
        groups: {
            ion: 0.0
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        }

    b3_ki_temp = {
        groups: {
            ion: 1.2 #b3_ki_full *NOTE:assumed to be constant in Zuend et al paper
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        }
    b3_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_abundance)),)
    b3_ki_temp_persistent[:,:]=1.2
    b1_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b1[main_group][ion-1]
            for ion, count in m_ions.items() #remember this 'ion' is numeric not a pybel object
            }
        for groups in non_zero_groups_old
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }
    b1_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_abundance)),)
    for groups in non_zero_groups_old:
        for ion_pybel, count in m_inorg.items():
            for ion, num in count.items():
                #pdb.set_trace()
                if num > 0:
                    b1_ki_temp_persistent[non_zero_groups_list.index(groups),ion_dict2array[ion_pybel]]=b1_ki_temp[groups][ion]

    #pdb.set_trace()
    b2_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b2[main_group][ion-1]
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }
    b2_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_abundance)),)
    #pdb.set_trace()
    for groups in non_zero_groups_old:
        for ion_pybel, count in m_inorg.items():
            for ion, num in count.items():
                if num > 0:
                    #pdb.set_trace()
                    b2_ki_temp_persistent[non_zero_groups_list.index(groups),ion_dict2array[ion_pybel]]=b2_ki_temp[groups][ion]

    persistant_data['b1_ki_temp_persistent']=b1_ki_temp_persistent
    persistant_data['b2_ki_temp_persistent']=b2_ki_temp_persistent
    persistant_data['b3_ki_temp_persistent']=b3_ki_temp_persistent

    #pdb.set_trace()
    
    ##b1ca, b2ca, b3ca
    b1_ca_temp = {
        cation: {
            anion: 0.0
            for anion in anions
            }
        for cation in cations
        }
    b2_ca_temp = {
        cation: {
            anion: 0.0
            for anion in anions
            }
        for cation in cations
        }
    b3_ca_temp = {
        cation: {
            anion: 0.0
            for anion in anions
            }
        for cation in cations
        }
    b1_ca_temp = {
        cation: {
            anion: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b1[cation-1][anion-1]
            for anion in anions
            }
        for cation in cations
        }
    b2_ca_temp = {
        cation: {
            anion: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b2[cation-1][anion-1]
            for anion in anions
            }
        for cation in cations
        }
    b3_ca_temp = {
        cation: {
            anion: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b3[cation-1][anion-1]
            for anion in anions
            }
        for cation in cations
        }
    #pdb.set_trace()
    b1_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    b2_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    b3_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    #pdb.set_trace()
    for cation in cations:
        for anion in anions:
            b1_ca_temp_persistent[cation2array[cation],anion2array[anion]]=b1_ca_temp[cation][anion]
            b2_ca_temp_persistent[cation2array[cation],anion2array[anion]]=b2_ca_temp[cation][anion]
            b3_ca_temp_persistent[cation2array[cation],anion2array[anion]]=b3_ca_temp[cation][anion]

    persistant_data['b1_ca_temp_persistent']=b1_ca_temp_persistent
    persistant_data['b2_ca_temp_persistent']=b2_ca_temp_persistent
    persistant_data['b3_ca_temp_persistent']=b3_ca_temp_persistent

    ##c1ca, c2ca, c3ca
    c1_ca_temp = {
        cation: {
            anion: 0.0 #data.AIOMFAC_MR_ION_ION_INTERACTIONS_c1[cations][anions]
            for anion in anions
            }
        for cation in cations
        }
    c1_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    c1_ca_temp_persistent[:,:]=0.0

    c2_ca_temp = {
        cation: {
            anion: 0.0  #populate default value and change in following if needs be
            for anion in anions
            }
        for cation in cations
        }

    c2_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    c2_ca_temp_persistent[:,:]=0.0
    c2_ca_temp = {
        cation: {
            anion: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c2[cation-1][anion-1]
            for anion in anions
            }
        for cation in cations
        }
    for cation in cations:
        for anion in anions:
            c2_ca_temp_persistent[cation2array[cation],anion2array[anion]]=c2_ca_temp[cation][anion]    

    c3_ca_temp_persistent=np.zeros((len(cations),len(anions)),)
    c3_ca_temp_persistent[:,:]=0.6
    c3_ca_temp = {
        cation: {
            anion: 0.6 #populate default value and change in following if needs be
            for anion in anions
            }
        for cation in cations
        }
    c3_ca_temp = {
        cation: {
            anion: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c3[cation-1][anion-1]
            for anion in anions
            }
        for cation in cations
        }
    for cation in cations:
        for anion in anions:
            c3_ca_temp_persistent[cation2array[cation],anion2array[anion]]=c3_ca_temp[cation][anion]    

    persistant_data['c1_ca_temp_persistent']=c1_ca_temp_persistent
    persistant_data['c2_ca_temp_persistent']=c2_ca_temp_persistent
    persistant_data['c3_ca_temp_persistent']=c3_ca_temp_persistent

    #Cation-cation interaction parameter matrix
    #for clarity, update all possible ions to have 0.0 value so this can be changed
    cation_dict={202:0.0,203:0.0,204:0.0,205:0.0,221:0.0,223:0.0}
    anion_dict={242:0.0,245:0.0,261:0.0}

    Rc_cdash = {}
    Rc_cdash = {
        cation: {
            cation_dash:value
            for cation_dash, value in cation_dict.items()
            }
        for cation, value in cation_dict.items()
        }
    Rc_cdash[204][205]=-0.154486414200000e+00
    Rc_cdash[205][204]=-0.154486414200000e+00
    #Persistent value
    Rc_cdash_persistent=np.zeros((len(cations),len(cations)),)
    for cation, value in cation_dict.items():
        for cation_dash, value in cation_dict.items():
            if cation in cation2array.keys() and cation_dash in cation2array.keys():
                Rc_cdash_persistent[cation2array[cation],cation2array[cation_dash]]=Rc_cdash[cation][cation_dash]
            
    persistant_data['Rc_cdash_persistent']=Rc_cdash_persistent

    #Molecular weight of seperate functional groups
    Mk_sol = {
        groups: data.AIOMFAC_MASS[groups]
        for groups in non_zero_groups_list
        }
    Mk_sol_persistent = np.zeros((len(non_zero_groups_list)),)
    for groups in non_zero_groups_list:
        Mk_sol_persistent[non_zero_groups_list.index(groups)]=Mk_sol[groups]
        
    persistant_data['Mk_sol_persistent']=Mk_sol_persistent

    Qc_cdash = {
        cation: {
            cation_dash:{
                anion: 0.0
                for anion, value in anion_dict.items()
                }
            for cation_dash, value in cation_dict.items()
            }
        for cation, value in cation_dict.items()
        }
    Qc_cdash[204][205][248]=0.448354085300000e-03
    Qc_cdash[205][204][248]=0.448354085300000e-03
    Qc_cdash_persistent=np.zeros((len(cations),len(cations),len(anions)),)
    for cation, value in cation_dict.items():
        for cation_dash, value in cation_dict.items():
            for anion, value in anion_dict.items():
                if cation in cation2array.keys() and cation_dash in cation2array.keys() and anion in anion2array.keys():
                    Qc_cdash_persistent[cation2array[cation],cation2array[cation_dash],anion2array[anion]]=Qc_cdash[cation][cation_dash][anion]
                    
    
    persistant_data['Qc_cdash_persistent']=Qc_cdash_persistent

    return persistant_data
    
def aiomfac_mr_quick(abundance_array,persistent_data,cation_index,anion_index,temperature):
    #This function relies on the provision of persistent variables and uses numpy operation wherever possible.
    
    
    non_zero_groups_list=persistent_data['non_zero_groups_list']
    non_zero_groups_flag=persistent_data['non_zero_groups_flag']
    num_orgs=persistent_data['num_orgs']
    molw_array_org=persistent_data['molw_array_org']
    cation2array=persistent_data['cation2array'] 
    anion2array=persistent_data['anion2array']
    ion_charge=persistent_data['ion_charge']
    b1_ki_temp_persistent=persistent_data['b1_ki_temp_persistent']
    b2_ki_temp_persistent=persistent_data['b2_ki_temp_persistent']
    b3_ki_temp_persistent=persistent_data['b3_ki_temp_persistent']
    b1_ca_temp_persistent=persistent_data['b1_ca_temp_persistent']
    b2_ca_temp_persistent=persistent_data['b2_ca_temp_persistent']
    b3_ca_temp_persistent=persistent_data['b3_ca_temp_persistent']   
    c1_ca_temp_persistent=persistent_data['c1_ca_temp_persistent']
    c2_ca_temp_persistent=persistent_data['c2_ca_temp_persistent']
    c3_ca_temp_persistent=persistent_data['c3_ca_temp_persistent']    
    Rc_cdash_persistent=persistent_data['Rc_cdash_persistent']
    Mk_sol_persistent=persistent_data['Mk_sol_persistent']    
    Qc_cdash_persistent=persistent_data['Qc_cdash_persistent']

    ion_abundance=abundance_array[num_orgs::]
    solvent_kg=np.sum(molw_array_org*np.array(abundance_array[0:num_orgs,0]))*1.0e-3
    #For the persistent version, here we just multiply the non_zero_groups_flag by the 
    #concentration of each compound
    conc_non_zero_groups = np.sum(np.multiply(non_zero_groups_flag,np.array(abundance_array[0:num_orgs])),axis=0)
    total_non_zero_groups = np.sum(conc_non_zero_groups)
    #Persistent version
    mole_frac_non_zero_groups = conc_non_zero_groups / total_non_zero_groups

    ion_molalities = ion_abundance / solvent_kg
    
    cation_molality = ion_abundance[[cation_index]]/solvent_kg
    anion_molality = ion_abundance[[anion_index]]/solvent_kg
    
    #pdb.set_trace()

    Ionic_strength = np.sum(np.multiply(np.power(ion_charge,2.0),ion_molalities))*0.5
    #print("solvent_kg",solvent_kg)
    #---------Now perform the interaction parameter calculations----------------
    #pdb.set_trace()
    #Organic-ion interactions
    if (Ionic_strength<250.0):
        
        temp=(np.exp(-1.0*c3_ca_temp_persistent*np.power(Ionic_strength,0.5)))
        temp2=np.multiply(c3_ca_temp_persistent,np.power(Ionic_strength,-0.5)*temp)
        Cca_dash_persistent=np.multiply(-0.5*c2_ca_temp_persistent,temp2)
        #pdb.set_trace()
        
        temp = np.exp(-1.0*b3_ki_temp_persistent*np.power(Ionic_strength,0.5))
        Bki_persistent = b1_ki_temp_persistent+ np.multiply(b2_ki_temp_persistent,np.exp(-1.0*b3_ki_temp_persistent*np.power(Ionic_strength,0.5)))
        Bki_dash_persistent = np.multiply(-0.5*b2_ki_temp_persistent,np.multiply(b3_ki_temp_persistent,np.power(Ionic_strength,-0.5)*temp))
        temp=np.exp(-1.0*b3_ca_temp_persistent*(np.power(Ionic_strength,0.5)))
        Bca_persistent=b1_ca_temp_persistent+np.multiply(b2_ca_temp_persistent,temp)
        
        #pdb.set_trace()

        temp=(np.exp(-1.0*b3_ca_temp_persistent*np.power(Ionic_strength,0.5)))
        temp2=np.multiply(b3_ca_temp_persistent,np.power(Ionic_strength,-0.5)*temp)
        #temp3=np.multiply(temp2,temp)
        #pdb.set_trace()
        Bca_dash_persistent=np.multiply(-0.5*b2_ca_temp_persistent,temp2)
        
        #pdb.set_trace()

        temp=np.exp(-1.0*c3_ca_temp_persistent*(np.power(Ionic_strength,0.5)))
        Cca_persistent=c1_ca_temp_persistent+np.multiply(c2_ca_temp_persistent,temp)
        
        #pdb.set_trace()

        temp=(np.exp(-1.0*c3_ca_temp_persistent*np.power(Ionic_strength,0.5)))
        temp2=np.multiply(c3_ca_temp_persistent,np.power(Ionic_strength,-0.5)*temp)
        Cca_dash_persistent=np.multiply(-0.5*c2_ca_temp_persistent,temp2)
        
        #pdb.set_trace()
        
    elif (Ionic_strength>=250.0):

        Bki_persistent = b1_ki_temp_persistent
        Bki_dash_persistent = np.zeros((len(non_zero_groups_list),len(ion_molalities)),)

        Bca_persistent=b1_ca_temp_persistent
        Bca_dash_persistent=np.zeros((len(cation_molality),len(anion_molality)),)

        Cca_persistent=c1_ca_temp_persistent
        Cca_dash_persistent=np.zeros((len(cation_molality),len(anion_molality)),)
    
    #pdb.set_trace()

    #Average molecular weight of the functional groups
    M_solv_mix_persistent = np.sum(Mk_sol_persistent*mole_frac_non_zero_groups)
    #Now calculate the activity coefficient of the separate organic functional groups
    summation1_persistent = np.sum(Bki_persistent*np.transpose(ion_molalities),axis=1)
    matrix_temp=np.repeat((mole_frac_non_zero_groups.reshape(len(non_zero_groups_list),1)),len(ion_molalities),axis=1)*np.transpose(ion_molalities)
    temp_summation1_persistent = np.sum(np.multiply((Bki_persistent+Ionic_strength*Bki_dash_persistent),matrix_temp))
    summation2_persistent=temp_summation1_persistent*(Mk_sol_persistent/M_solv_mix_persistent)
    
    #pdb.set_trace()

    #A2) ion-ion interactions
    matrix_temp = np.repeat((cation_molality.reshape(len(cation_molality),1)),len(anion_molality),axis=1)*anion_molality
    temp_summation2_persistent = np.sum(np.multiply((Bca_persistent + Ionic_strength*Bca_dash_persistent),matrix_temp))
    summation3_persistent=temp_summation2_persistent*(Mk_sol_persistent)

    #matrix_temp = np.repeat((cation_molality.reshape(len(cation_molality),1)),len(anion_molality),axis=1)*anion_molality
    #temp_summation2_persistent = np.sum(np.multiply((Bca_persistent + Ionic_strength*Bca_dash_persistent),matrix_temp))
    #summation3_persistent=temp_summation2_persistent*(Mk_sol_persistent)

    #A3) ion-ion interactions
    temp_summation3_persistent = np.sum(np.multiply(ion_molalities,abs(ion_charge)))
    temp_summation4_persistent = np.sum(np.multiply((2.0*Cca_persistent + Ionic_strength*Cca_dash_persistent),matrix_temp))
    summation4_persistent=temp_summation3_persistent*temp_summation4_persistent*(Mk_sol_persistent)
    
    #pdb.set_trace()
    
    #A4) cation, cation (only need to do if more than one cation)
    matrix_temp2 = np.repeat((cation_molality.reshape(len(cation_molality),1)),len(cation_molality),axis=1)*cation_molality
    temp_summation5_persistent = np.sum(np.multiply(Rc_cdash_persistent,matrix_temp2))
    summation5_persistent = temp_summation5_persistent*Mk_sol_persistent
    temp_summation6_persistent = 0.0
    step=0
    for value in anion_molality:
        temp_summation6_persistent+=value*np.sum(np.multiply(2.0*Qc_cdash_persistent[:,:,step],matrix_temp2))
        step+=1
    summation6_persistent = temp_summation6_persistent*Mk_sol_persistent

    ln_gamma_k_MR_persistent = summation1_persistent - summation2_persistent - summation3_persistent - summation4_persistent - summation5_persistent - summation6_persistent
    #pdb.set_trace()
    #pdb.set_trace()
    #print(ln_gamma_k_MR_persistent)
    #Ln_gamma_s_MR=numpy.zeros((1,org_molecules),)
    #for molecule_step in range(org_molecules):
    #    for k in range(max_group_num_org_main):
    #        Ln_gamma_s_MR[0,molecule_step]=Ln_gamma_s_MR[0,molecule_step]+org_group_stoich_new_main[molecule_step,k]*ln_gamma_k_MR[0,k]

    #Now sum all of the individual functional groups together to get the final activity coefficient
    #Add up contributions from each group according to stochiometry
    Ln_gamma_s_MR_persistent =np.sum((non_zero_groups_flag*ln_gamma_k_MR_persistent),axis=1)
    
    #pdb.set_trace()
    
    #Activity coefficient for ions
    #Generic calculations
    summation1_ion_persistent = np.sum(Bki_persistent*mole_frac_non_zero_groups[:, np.newaxis],axis=0)*(1.0/M_solv_mix_persistent)
    
    #pdb.set_trace()

    pre_summation2_ion_persistent = np.sum(np.multiply(np.multiply(Bki_dash_persistent,np.transpose(ion_molalities)),mole_frac_non_zero_groups[:, np.newaxis]))
    #print("pre_summation2_ion_persistent [python]",pre_summation2_ion_persistent)
    #print("M_solv_mix_persistent [python]",M_solv_mix_persistent)

    summation2_ion_persistent = (np.power(abs(ion_charge),2.0)/2.0*M_solv_mix_persistent)*pre_summation2_ion_persistent
    
    #pdb.set_trace()

    pre_summation4_ion_persistent = np.sum(matrix_temp*Bca_dash_persistent)
    summation4_ion_persistent = pre_summation4_ion_persistent*0.5*np.power(abs((ion_charge)),2.0)
    
    #pdb.set_trace()

    pre_summation5_ion_persistent = np.sum(abs(ion_charge)*ion_molalities)
    pre_summation6_ion_persistent = pre_summation5_ion_persistent
    summation6_ion_persistent=np.zeros((len(ion_molalities)),)
    step=0
    for value in ion_molalities:
        summation6_ion_persistent[step]=np.sum(((Cca_persistent*abs(ion_charge[step]))+Cca_dash_persistent*(np.power(abs(ion_charge[step]),2.0)*0.5)*pre_summation6_ion_persistent)*matrix_temp)
        step+=1
    
    #pdb.set_trace()
    
    #Cation/anion specific calculations
    summation3_ion_cation_persistent = np.sum(np.multiply(Bca_persistent,anion_molality[np.newaxis,:]),axis=1)
    summation3_ion_anion_persistent = np.sum(np.multiply(Bca_persistent,cation_molality[:,np.newaxis]),axis=0)
    
    #pdb.set_trace()

    summation3_ion_persistent =summation3_ion_cation_persistent.copy()
    summation3_ion_persistent=np.append(summation3_ion_persistent,summation3_ion_anion_persistent)
    
    #pdb.set_trace()

    summation5_ion_cation_persistent=np.sum(np.multiply(Cca_persistent,anion_molality[np.newaxis,:])*pre_summation5_ion_persistent,axis=1)
    summation5_ion_anion_persistent=np.sum(np.multiply(Cca_persistent,cation_molality[:,np.newaxis])*pre_summation5_ion_persistent,axis=0)
    
    #pdb.set_trace()

    summation5_ion_persistent =summation5_ion_cation_persistent.copy()
    summation5_ion_persistent=np.append(summation5_ion_persistent,summation5_ion_anion_persistent)
    
    #pdb.set_trace()
    
    summation7_ion_cation_1_persistent =np.sum(np.multiply(Rc_cdash_persistent,cation_molality[np.newaxis,:]),axis=1)
    
    #pdb.set_trace()

    summation7_ion_cation_2_persistent=np.zeros((len(cation_molality)),)
    step=0
    for value in cation_molality:
        summation7_ion_cation_2_persistent[step]=np.sum(np.multiply(Qc_cdash_persistent[step,:,:],matrix_temp))
        
    #pdb.set_trace()
        
    summation7_ion_cation_persistent=summation7_ion_cation_1_persistent+summation7_ion_cation_2_persistent
    
    #pdb.set_trace()

    summation7_ion_anion_persistent=np.zeros((len(anion_molality)),)
    step=0
    for value in anion_molality:
        summation7_ion_anion_persistent[step]=np.sum(np.multiply(Qc_cdash_persistent[:,:,step],matrix_temp))
        
    #pdb.set_trace()

    summation7_ion_persistent=summation7_ion_cation_persistent.copy()
    summation7_ion_persistent=np.append(summation7_ion_persistent,summation7_ion_anion_persistent)
    
    #print("summation1_ion_persistent [python]",summation1_ion_persistent)
    #print("summation2_ion_persistent [python]",summation2_ion_persistent)
    #print("summation3_ion_persistent [python]",summation3_ion_persistent)
    #print("summation4_ion_persistent [python]",summation4_ion_persistent)
    #print("summation5_ion_persistent [python]",summation5_ion_persistent)
    #print("summation6_ion_persistent [python]",summation6_ion_persistent)
    #print("summation7_ion_persistent [python]",summation7_ion_persistent)

    #pdb.set_trace()
    
    Ln_gamma_i_MR_persistent=summation1_ion_persistent+summation2_ion_persistent[:,0]+summation3_ion_persistent+summation4_ion_persistent[:,0]+\
    summation5_ion_persistent+summation6_ion_persistent+summation7_ion_persistent
    
    #pdb.set_trace()

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

    #pdb.set_trace()

    #--Organic solvent activity coefficients---
    Ln_gamma_s_LR_persistent=((2.0*A*(molw_array_org*1e-3))/(np.power(b,3.0)))*summation_org
    
    #pdb.set_trace()

    #--Ionic activity coefficients--
    Ln_gamma_i_LR_persistent=(-1.0*np.power(abs(ion_charge),2.0)*summation_ion_1)/summation_ion_2
    #print("summation_ion_2 [python]",summation_ion_2)

    #print("Ln_gamma_i_LR_persistent = ", Ln_gamma_i_LR_persistent)
    #pdb.set_trace()

    #Now calculate the conversion factor to change the reference state from mole fraction to molality    
    x_i_org_persistent=np.array(abundance_array[0:num_orgs,0])/(np.sum(np.array(abundance_array[0:num_orgs,0])))
    
    #pdb.set_trace()

    mean_mw_solvent_persistent=np.sum(np.multiply(x_i_org_persistent,molw_array_org)*1e-3)
    
    #pdb.set_trace()

    sum_molalities_persistent=np.sum(ion_molalities)
    
    #pdb.set_trace()

    ionic_conversion_factor_persistent=log((0.01801528/mean_mw_solvent_persistent)+0.01801528*sum_molalities_persistent)
    #print("ionic_conversion_factor_persistent [python] = ", ionic_conversion_factor_persistent,)
    #pdb.set_trace()

    #Now map this ontp pybel objects as keys rather than integers
    Ln_gamma_i_MR_LR_persistent=Ln_gamma_i_MR_persistent+Ln_gamma_i_LR_persistent[:,0]-ionic_conversion_factor_persistent
    
    #pdb.set_trace()

    #Ln_gamma_i_MR_LR_keys = {
    #    ion: Ln_gamma_i_MR_LR[old_key] #q_k_i[compound].get(group1, 0)
    #    for ion, matches in m_inorg.items()
    #    for old_key, count in matches.items()
    #    if count >0
    #    }
    #Ln_gamma_s_MR_LR = {
    #    compound: Ln_gamma_s_LR[compound]+Ln_gamma_s_MR[compound]
    #    for compound, abundance in organic_compounds.items()
    #    }
    Ln_gamma_s_MR_LR_persistent = Ln_gamma_s_LR_persistent+Ln_gamma_s_MR_persistent
    Ln_gamma_tot_MR_LR_persistent = np.append(Ln_gamma_s_MR_LR_persistent.T,Ln_gamma_i_MR_LR_persistent)
    
    #pdb.set_trace()
        
    # Create a dictionary that holds both ionic and organic values
    #Ln_gamma_tot_MR_LR={}
    #for compound, abundance in organic_compounds.items():
    #    Ln_gamma_tot_MR_LR[compound] = Ln_gamma_s_MR_LR[compound]
    # 
    #for ion, matches in m_inorg.items():
    #    Ln_gamma_tot_MR_LR[ion]= Ln_gamma_i_MR_LR_keys[ion]


    #NOTE we only return the organic activity coefficients here. For including inorganic
    #ions, this will also need to be passed
    #pdb.set_trace()
    return Ln_gamma_tot_MR_LR_persistent 
    
    
def aiomfac_mr_test(species_dict2array,ion_dict2array,num_species,abundance_array,cation_index,anion_index,organic_compounds, inorganic_ions, temperature):

    m = aiomfac_organic(organic_compounds)
    m.update(aiomfac_inorganic(inorganic_ions))

    m_org = aiomfac_organic(organic_compounds)
    m_inorg = aiomfac_inorganic(inorganic_ions)

    list_compounds=organic_compounds.copy()
    list_compounds.update(inorganic_ions)
    #Mole fraction of all compounds

    #For the persistent version we need some way of extracting ions from the abundance array
    #We will have to do this from the calling script.
    #pdb.set_trace()
    ion_abundance=abundance_array[num_species::]

    non_zero_groups_old = {
        group
        for compound, matches in m_org.items()
        for group, count in matches.items()
        if count > 0
        }
    #--Persistent variable--
    non_zero_groups_list = []
    for group in non_zero_groups_old:
        if group not in non_zero_groups_list:
            non_zero_groups_list.append(group)

    #pdb.set_trace()

    #Now create a matrix that simply records the flag on a nonzero group or not
    #Persistent variable
    non_zero_groups_flag_old = {
        compound : { group : count
        for group, count in matches.items()
        }
        for compound, matches in m_org.items()
    }
    #pdb.set_trace()
    #--Persistent variable--
    non_zero_groups_flag = np.zeros((len(m_org.keys()),len(non_zero_groups_list)),)
    for compound, matches in m_org.items():
        for group, count in matches.items():
            if group in non_zero_groups_list:
                non_zero_groups_flag[species_dict2array[compound],non_zero_groups_list.index(group)]=count

    #pdb.set_trace()
    ##-------------TESTING AIOMFAC MR TERMS------------------

    solvent_kg_old=sum((key.molwt*1e-3)*value for key, value in organic_compounds.items())
    #--Persistent variable--
    #Persistent variable - [molw_array]
    num_orgs=len(organic_compounds.keys())
    #--Persistent variable--
    molw_array_org=np.zeros((num_orgs),)
    #pdb.set_trace()
    for key, value in organic_compounds.items():
        molw_array_org[species_dict2array[key]]=key.molwt
    #Non persistent variable
    #pdb.set_trace()
    solvent_kg=np.sum(molw_array_org*np.array(abundance_array[0:num_orgs,0]))*1.0e-3

    #pdb.set_trace()

    #total_non_zero_groups=0
    #conc_non_zero_groups={}
    #for compound, matches in m_org.items():
    #    for group, count in matches.items():
    #        if count > 0:
    #            #conc_non_zero_groups[group]=conc_non_zero_groups[group]+count
    #            total_non_zero_groups=total_non_zero_groups+count

    #only over the organic functional groups
    conc_non_zero_groups_old = {
        group: sum(m_org[compound].get(group, 0)*organic_compounds[compound] for compound, matches in m_org.items())
        for group in non_zero_groups_old
                           }
    
    #--Persistent variable--
    #For the persistent version, here we just multiply the non_zero_groups_flag by the 
    #concentration of each compound
    conc_non_zero_groups = np.sum(np.multiply(non_zero_groups_flag,np.array(abundance_array[0:num_orgs])),axis=0)
    total_non_zero_groups = np.sum(conc_non_zero_groups)
    
                           
    total_non_zero_groups_old=sum(conc_non_zero_groups_old.values())                       
    #again, only over the organic compounds
    mole_frac_non_zero_groups_old = {
        group: conc_non_zero_groups_old[group]/total_non_zero_groups_old
        for group in non_zero_groups_old
                                }
    #pdb.set_trace()
    #Persistent version
    mole_frac_non_zero_groups = conc_non_zero_groups / total_non_zero_groups
    #pdb.set_trace()
    #now calculate properties associated with the ions - molalities
    m_ions = groups_func.all_matches(data.AIOMFAC_ION_SMARTS, inorganic_ions) #might use this
  
    #pdb.set_trace()
    ion_molalities_old = {ions: value/solvent_kg_old for ions, value in m_ions.iteritems()}
    #pdb.set_trace()
    #For the persistent version transfer the keys to matrix entries 
    #the ions abundance will be tramsferred outside of this loop
    #Persistent version - ions ordered according to fixed addition
    #We also pass in :
    #cation/anion index from the initialisation of the model
    #ion charge array
    #pdb.set_trace()
    ion_molalities = ion_abundance / solvent_kg

    #Cation/Anion molality
    cations = {
        ion
        for ion, count in m_ions.items() #ion is a numeric key, not a pybel object
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    #Now setup a dictionary that creates cation array indices
    step=0
    cation2array = dict()
    for ion in cations:
        cation2array[ion]=step
        step+=1

    cation_molality_old = {
        ion: m_ions[ion]/solvent_kg_old
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] > 0.0
        }
    
    cation_molality = ion_abundance[[cation_index]]/solvent_kg

    anions  = {
        ion
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    step=0
    anion2array = dict()
    for ion in anions:
        anion2array[ion]=step
        step+=1

    anion_molality_old = {
        ion: m_ions[ion]/solvent_kg_old
        for ion, count in m_ions.items()
        if count and data.AIOMFAC_ION_CHARGE[ion] < 0.0
        }
    anion_molality = ion_abundance[[anion_index]]/solvent_kg

    #Ionic strength
    Ionic_strength = 0
    for c, value in cation_molality_old.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[c], 2)))
    for a, value in anion_molality_old.items():
        Ionic_strength = Ionic_strength + (value*(pow(data.AIOMFAC_ION_CHARGE[a], 2)))
    Ionic_strength = Ionic_strength*0.5

    #Persistent variable
    ion_charge=np.zeros((len(ion_abundance),1),)
    for ion, info in m_inorg.items():
        for key, value in info.iteritems():
            if value > 0:
                key_rec=key
                ion_charge[ion_dict2array[ion],0]=data.AIOMFAC_ION_CHARGE[key_rec]
    #pdb.set_trace()    
    #Persistent version
    Ionic_strength_new = np.sum(np.multiply(np.power(ion_charge,2.0),ion_molalities))*0.5

    
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
        for groups in non_zero_groups_old
        }

    b2_ki_temp = {
        groups: {
            ion: 0.0
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        }

    b3_ki_temp = {
        groups: {
            ion: 1.2 #b3_ki_full *NOTE:assumed to be constant in Zuend et al paper
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        }
    b3_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_molalities)),)
    b3_ki_temp_persistent[:,:]=1.2
    b1_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b1[main_group][ion-1]
            for ion, count in m_ions.items() #remember this 'ion' is numeric not a pybel object
            }
        for groups in non_zero_groups_old
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }
    b1_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_molalities)),)
    for groups in non_zero_groups_old:
        for ion_pybel, count in m_inorg.items():
            for ion, num in count.iteritems():
                #pdb.set_trace()
                if num > 0:
                    b1_ki_temp_persistent[non_zero_groups_list.index(groups),ion_dict2array[ion_pybel]]=b1_ki_temp[groups][ion]

    #pdb.set_trace()
    b2_ki_temp = {
        groups: {
            ion: data.AIOMFAC_MR_ORG_ION_INTERACTIONS_b2[main_group][ion-1]
            for ion, count in m_ions.items()
            }
        for groups in non_zero_groups_old
        for main_group in (data.AIOMFAC_MAIN_GROUP[groups],)
        }
    b2_ki_temp_persistent = np.zeros((len(non_zero_groups_list),len(ion_molalities)),)
    #pdb.set_trace()
    for groups in non_zero_groups_old:
        for ion_pybel, count in m_inorg.items():
            for ion, num in count.iteritems():
                if num > 0:
                    #pdb.set_trace()
                    b2_ki_temp_persistent[non_zero_groups_list.index(groups),ion_dict2array[ion_pybel]]=b2_ki_temp[groups][ion]

    
    ##b1ca, b2ca, b3ca
    b1_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    b2_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    b3_ca_temp = {
        cations: {
            anions: 0.0
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    b1_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b1[cations-1][anions-1]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    b2_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b2[cations-1][anions-1]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    b3_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_b3[cations-1][anions-1]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    #pdb.set_trace()
    b1_ca_temp_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)
    b2_ca_temp_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)
    b3_ca_temp_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)
    #pdb.set_trace()
    for cations, value in cation_molality_old.items():
        for anions, value2 in anion_molality_old.items():
            b1_ca_temp_persistent[cation2array[cations],anion2array[anions]]=b1_ca_temp[cations][anions]
            b2_ca_temp_persistent[cation2array[cations],anion2array[anions]]=b2_ca_temp[cations][anions]
            b3_ca_temp_persistent[cation2array[cations],anion2array[anions]]=b3_ca_temp[cations][anions]

    
    ##c1ca, c2ca, c3ca
    c1_ca_temp = {
        cations: {
            anions: 0.0 #data.AIOMFAC_MR_ION_ION_INTERACTIONS_c1[cations][anions]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    c1_ca_temp_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)
    c1_ca_temp_persistent[:,:]=0.0

    c2_ca_temp = {
        cations: {
            anions: 0.0  #populate default value and change in following if needs be
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }

    c2_ca_temp_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)
    c2_ca_temp_persistent[:,:]=0.0
    c2_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c2[cations-1][anions-1]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    for cations, value in cation_molality_old.items():
        for anions, value in anion_molality_old.items():
            c2_ca_temp_persistent[cation2array[cations],anion2array[anions]]=c2_ca_temp[cations][anions]    

    c3_ca_temp_persistent=np.zeros((len(cation_molality_old.items()),len(anion_molality_old.items())),)
    c3_ca_temp_persistent[:,:]=0.6
    c3_ca_temp = {
        cations: {
            anions: 0.6 #populate default value and change in following if needs be
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    c3_ca_temp = {
        cations: {
            anions: data.AIOMFAC_MR_ION_ION_INTERACTIONS_c3[cations-1][anions-1]
            for anions, value in anion_molality_old.items()
            }
        for cations, value in cation_molality_old.items()
        }
    for cations, value in cation_molality_old.items():
        for anions, value in anion_molality_old.items():
            c3_ca_temp_persistent[cation2array[cations],anion2array[anions]]=c3_ca_temp[cations][anions]    


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
    #Persistent value
    Rc_cdash_persistent=np.zeros((len(cation_molality_old.keys()),len(cation_molality_old.keys())),)
    for cations, value in cation_dict.items():
        for cations_dash, value in cation_dict.items():
            if cations in cation2array.keys() and cations_dash in cation2array.keys():
                Rc_cdash_persistent[cation2array[cations],cation2array[cations_dash]]=Rc_cdash[cations][cations_dash]

    #Molecular weight of seperate functional groups
    Mk_sol = {
        groups: data.AIOMFAC_MASS[groups]
        for groups in non_zero_groups_list
        }
    Mk_sol_persistent = np.zeros((len(non_zero_groups_list)),)
    for groups in non_zero_groups_list:
        Mk_sol_persistent[non_zero_groups_list.index(groups)]=Mk_sol[groups]

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
    Qc_cdash_persistent=np.zeros((len(cation_molality_old.items()),len(cation_molality_old.items()),len(anion_molality_old.items())),)
    for cations, value in cation_molality_old.items():
        for cations_dash, value in cation_molality_old.items():
            for anions, value in anion_molality_old.items():
                if cations in cation2array.keys() and cations_dash in cation2array.keys() and anions in anion2array.keys():
                    Qc_cdash_persistent[cation2array[cations],cation2array[cations_dash],anion2array[anions]]=Qc_cdash[cations][cations_dash][anions]

    # Need to perform code check up to here!!
    
    #---------Now perform the interaction parameter calculations----------------
    #pdb.set_trace()
    #Organic-ion interactions
    if (Ionic_strength<250.0):
        Bki = {
            groups: {
                ions : b1_ki_temp[groups][ions]+b2_ki_temp[groups][ions]*exp(-1.0*b3_ki_temp[groups][ions]*(pow(Ionic_strength,0.5)))
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups_list
            }
        #pdb.set_trace()
        Bki_persistent = b1_ki_temp_persistent+ np.multiply(b2_ki_temp_persistent,np.exp(-1.0*b3_ki_temp_persistent*np.power(Ionic_strength,0.5)))
        Bki_dash = {
            groups: {
                #Bki_dash[k,i]=-0.5*b2_ki_temp[k,i]*b3_ki_temp[k,i]*(numpy.power(I,-0.5))*numpy.exp(-1*b3_ki_temp[k,i]*numpy.power(I,0.5))
                ions: -0.5*b2_ki_temp[groups][ions]*b3_ki_temp[groups][ions]*(pow(Ionic_strength,-0.5))*(exp(-1*b3_ki_temp[groups][ions]*pow(Ionic_strength,0.5)))
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups_list
            }
        temp = np.exp(-1.0*b3_ki_temp_persistent*np.power(Ionic_strength,0.5))
        #pdb.set_trace()
        Bki_dash_persistent = np.multiply(-0.5*b2_ki_temp_persistent,np.multiply(b3_ki_temp_persistent,np.power(Ionic_strength,-0.5)*temp))
        Bca = {
            cations: {
                anions : b1_ca_temp[cations][anions]+b2_ca_temp[cations][anions]*exp(-1.0*b3_ca_temp[cations][anions]*(pow(Ionic_strength,0.5)))
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        temp=np.exp(-1.0*b3_ca_temp_persistent*(np.power(Ionic_strength,0.5)))
        #pdb.set_trace()
        Bca_persistent=b1_ca_temp_persistent+np.multiply(b2_ca_temp_persistent,temp)
        pdb.set_trace()
        Bca_dash = {
            cations: {
                #Bca_dash[c,a]=-0.5*b2_ca_temp[c,a]*b3_ca_temp[c,a]*(numpy.power(I,-0.5))*numpy.exp(-1*b3_ca_temp[c,a]*numpy.power(I,0.5))
                anions: -0.5*b2_ca_temp[cations][anions]*b3_ca_temp[cations][anions]*(pow(Ionic_strength,-0.5))*(exp(-1*b3_ca_temp[cations][anions]*pow(Ionic_strength,0.5)))
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        temp=(np.exp(-1.0*b3_ca_temp_persistent*np.power(Ionic_strength,0.5)))
        temp2=np.multiply(b3_ca_temp_persistent,np.power(Ionic_strength,-0.5)*temp)
        #temp3=np.multiply(temp2,temp)
        #pdb.set_trace()
        Bca_dash_persistent=np.multiply(-0.5*b2_ca_temp_persistent,temp2)
        Cca = {
            cations: {
                anions : c1_ca_temp[cations][anions]+c2_ca_temp[cations][anions]*exp(-1.0*c3_ca_temp[cations][anions]*(pow(Ionic_strength,0.5)))
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        temp=np.exp(-1.0*c3_ca_temp_persistent*(np.power(Ionic_strength,0.5)))
        #pdb.set_trace()
        Cca_persistent=c1_ca_temp_persistent+np.multiply(c2_ca_temp_persistent,temp)
        Cca_dash = {
            cations: {
                #Cca_dash[c,a]=-0.5*c2_ca_temp[c,a]*c3_ca_temp[c,a]*(numpy.power(I,-0.5))*numpy.exp(-1*c3_ca_temp[c,a]*numpy.power(I,0.5))
                anions: -0.5*c2_ca_temp[cations][anions]*c3_ca_temp[cations][anions]*(pow(Ionic_strength,-0.5))*(exp(-1*c3_ca_temp[cations][anions]*pow(Ionic_strength,0.5)))
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        temp=(np.exp(-1.0*c3_ca_temp_persistent*np.power(Ionic_strength,0.5)))
        temp2=np.multiply(c3_ca_temp_persistent,np.power(Ionic_strength,-0.5)*temp)
        Cca_dash_persistent=np.multiply(-0.5*c2_ca_temp_persistent,temp2)
        #pdb.set_trace()

    elif (Ionic_strength>=250.0):
        Bki = {
            groups: {
                ions : b1_ki_temp[groups][ions]
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups_list
            }
        Bki_persistent = b1_ki_temp_persistent
        Bki_dash = {
            groups: {
                ions: 0.0
                for ions, count in m_ions.items()
                }
            for groups in non_zero_groups_list
            }
        Bki_dash_persistent = np.zeros((len(non_zero_groups_list),len(ion_molalities)),)
        Bca = {
            cations: {
                anions : b1_ca_temp[cations][anions]
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        Bca_persistent=b1_ca_temp_persistent
        Bca_dash = {
            cations: {
                anions: 0.0
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        Bca_dash_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)

        Cca = {
            cations: {
                anions : c1_ca_temp[cations][anions]
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        Cca_persistent=c1_ca_temp_persistent
        Cca_dash = {
            cations: {
                anions: 0.0
                for anions, value in anion_molality_old.items()
                }
            for cations, value in cation_molality_old.items()
            }
        Cca_dash_persistent=np.zeros((len(cation_molality_old.keys()),len(anion_molality_old.keys())),)

    pdb.set_trace()
    #Average molecular weight of the functional groups
    M_solv_mix=sum(Mk_sol[key]*mole_frac_non_zero_groups_old[key] for key in non_zero_groups_list)
    pdb.set_trace()
    M_solv_mix_persistent = np.sum(Mk_sol_persistent*mole_frac_non_zero_groups)

    #Now calculate the activity coefficient of the separate organic functional groups
    summation1 = {
        groups: sum (Bki[groups][ions]*value for ions, value in ion_molalities_old.iteritems())
        for groups in non_zero_groups_list
        } 
    pdb.set_trace()
    summation1_persistent = np.sum(Bki_persistent*np.transpose(ion_molalities),axis=1)
    temp_summation1=sum((Bki[groups][ions]+Ionic_strength*Bki_dash[groups][ions])*\
            (mole_frac_non_zero_groups_old[groups]*value)
            for ions, value in ion_molalities_old.iteritems()
            for groups in non_zero_groups_list)
    #For the persistent value, we work across each row first, and then sum the total
    #We ned to first create a new matrix of mole_frac * molalities 
    #matrix_temp=np.zeros((len(non_zero_groups_list),len(ion_molalities)),)
    pdb.set_trace()
    matrix_temp=np.repeat((mole_frac_non_zero_groups.reshape(len(non_zero_groups_list),1)),len(ion_molalities),axis=1)*np.transpose(ion_molalities)
    temp_summation1_persistent = np.sum(np.multiply((Bki_persistent+Ionic_strength_new*Bki_dash_persistent),matrix_temp))
    summation2 = {
        groups:temp_summation1*(Mk_sol[groups]/M_solv_mix)
        for groups in non_zero_groups_list
        }
    pdb.set_trace()
    summation2_persistent=temp_summation1_persistent*(Mk_sol_persistent/M_solv_mix_persistent)

    #A2) ion-ion interactions
    temp_summation2=sum((Bca[cations][anions]+Ionic_strength*Bca_dash[cations][anions])*cation_molality_old[cations]*anion_molality_old[anions]
        for anions, value in anion_molality_old.iteritems()
        for cations, value in cation_molality_old.iteritems())

    matrix_temp = np.repeat((cation_molality.reshape(len(cation_molality),1)),len(anion_molality),axis=1)*anion_molality
    temp_summation2_persistent = np.sum(np.multiply((Bca_persistent + Ionic_strength*Bca_dash_persistent),matrix_temp))
    summation3 = {
        groups:temp_summation2*Mk_sol[groups]
        for groups in non_zero_groups_list
        }
    
    summation3_persistent=temp_summation2_persistent*(Mk_sol_persistent)
    #A3) ion-ion interactions
    temp_summation3 = 0
    for c, value in cation_molality_old.items():
        temp_summation3 = temp_summation3 + (value*(abs(data.AIOMFAC_ION_CHARGE[c])))
    for a, value in anion_molality_old.items():
        temp_summation3 = temp_summation3 + (value*(abs(data.AIOMFAC_ION_CHARGE[a])))
    pdb.set_trace()
    temp_summation3_persistent = np.sum(np.multiply(ion_molalities,abs(ion_charge)))
    temp_summation4 = sum ((2.0*Cca[cations][anions]+Ionic_strength*Cca_dash[cations][anions])*ion_molalities_old[cations]*ion_molalities_old[anions]
        for anions, value in anion_molality_old.iteritems()
        for cations, value in cation_molality_old.iteritems())
    pdb.set_trace()
    temp_summation4_persistent = np.sum(np.multiply((2.0*Cca_persistent + Ionic_strength*Cca_dash_persistent),matrix_temp))
    summation4 = {
        groups:temp_summation3*temp_summation4*Mk_sol[groups]
        for groups in non_zero_groups_list
        }
    pdb.set_trace()
    summation4_persistent=temp_summation3_persistent*temp_summation4_persistent*(Mk_sol_persistent)
    #A4) cation, cation (only need to do if more than one cation)
    temp_summation5 = sum (Rc_cdash[cations][cations_dash]*ion_molalities_old[cations]*ion_molalities_old[cations_dash]
        for cations_dash, value in cation_molality_old.items()
        for cations, value in cation_molality_old.items())
    matrix_temp2 = np.repeat((cation_molality.reshape(len(cation_molality),1)),len(cation_molality),axis=1)*cation_molality
    pdb.set_trace()
    temp_summation5_persistent = np.sum(np.multiply(Rc_cdash_persistent,matrix_temp2))
    summation5 = {
        groups:temp_summation5*Mk_sol[groups]
        for groups in non_zero_groups_list
        }
    summation5_persistent = temp_summation5_persistent*Mk_sol_persistent
    temp_summation6 = sum(2.0*Qc_cdash[cations][cations_dash][anions]*c_value*cdash_value*a_value
        for anions, a_value in anion_molality_old.items()
        for cations_dash, cdash_value in cation_molality_old.items()
        for cations, c_value in cation_molality_old.items())
    temp_summation6_persistent = 0.0
    step=0
    pdb.set_trace()
    for value in anion_molality:
        temp_summation6_persistent+=value*np.sum(np.multiply(2.0*Qc_cdash_persistent[:,:,step],matrix_temp2))
        step+=1
    summation6 = {
        groups:temp_summation6*Mk_sol[groups]
        for groups in non_zero_groups_list
        }
    summation6_persistent = temp_summation6_persistent*Mk_sol_persistent
    pdb.set_trace()
    ln_gamma_k_MR = {
        groups: summation1[groups]-summation2[groups]-summation3[groups]-summation4[groups]-summation5[groups]-summation6[groups]
        for groups in non_zero_groups_list
        }
    ln_gamma_k_MR_persistent = summation1_persistent - summation2_persistent - summation3_persistent - summation4_persistent - summation5_persistent - summation6_persistent
    #pdb.set_trace()
    #Ln_gamma_s_MR=numpy.zeros((1,org_molecules),)
    #for molecule_step in range(org_molecules):
    #    for k in range(max_group_num_org_main):
    #        Ln_gamma_s_MR[0,molecule_step]=Ln_gamma_s_MR[0,molecule_step]+org_group_stoich_new_main[molecule_step,k]*ln_gamma_k_MR[0,k]

    #Now sum all of the individual functional groups together to get the final activity coefficient
    pdb.set_trace()
    Ln_gamma_s_MR = {
        compound: sum(ln_gamma_k_MR[group]*m_org[compound][group] for group in non_zero_groups_list)
    for compound, abundance in organic_compounds.items()
    }
    pdb.set_trace()
    #Add up contributions from each group according to stochiometry
    Ln_gamma_s_MR_persistent =np.sum((non_zero_groups_flag*ln_gamma_k_MR_persistent),axis=1)
    pdb.set_trace()

    #Activity coefficient for ions
    #Generic calculations
    summation1_ion = {
    ion: sum(Bki[group][ion]*mole_frac_non_zero_groups_old[group]*(1.0/M_solv_mix)
    for group in non_zero_groups_list)
        for ion, value in ion_molalities_old.items()
    }
    pdb.set_trace()
    summation1_ion_persistent = np.sum(Bki_persistent*mole_frac_non_zero_groups[:, np.newaxis],axis=0)*(1.0/M_solv_mix_persistent)

    pre_summation2_ion = sum(Bki_dash[group][ion]*mole_frac_non_zero_groups_old[group]*value
    for ion, value in ion_molalities_old.items()
    for group in non_zero_groups_list)
    pdb.set_trace()
    pre_summation2_ion_persistent = np.sum(np.multiply(np.multiply(Bki_dash_persistent,np.transpose(ion_molalities)),mole_frac_non_zero_groups[:, np.newaxis]))

    summation2_ion = {
    ion:pre_summation2_ion*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)/2.0*M_solv_mix)
    for ion, value in ion_molalities_old.items()
    }
    pdb.set_trace()
    summation2_ion_persistent = (np.power(abs(ion_charge),2.0)/2.0*M_solv_mix_persistent)*pre_summation2_ion_persistent

    pre_summation4_ion = sum(
    Bca_dash[cation][anion]*c_val*a_val
    for anion, a_val in anion_molality_old.items()
    for cation, c_val in cation_molality_old.items()
    )
    pdb.set_trace()
    pre_summation4_ion_persistent = np.sum(matrix_temp*Bca_dash_persistent)
 
    summation4_ion = {
    ion: pre_summation4_ion*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)*0.5)
    for ion, value in ion_molalities_old.items()
    }
    pdb.set_trace()
    summation4_ion_persistent = pre_summation4_ion_persistent*0.5*np.power(abs((ion_charge)),2.0)
    pdb.set_trace()
    pre_summation5_ion=sum(abs(data.AIOMFAC_ION_CHARGE[ion])*value for ion, value in ion_molalities_old.items())
    pre_summation5_ion_persistent = np.sum(abs(ion_charge)*ion_molalities)
    pre_summation6_ion = pre_summation5_ion
    pre_summation6_ion_persistent = pre_summation5_ion_persistent

    summation6_ion = {
        ion:sum((Cca[cation][anion]*abs(data.AIOMFAC_ION_CHARGE[ion])+Cca_dash[cation][anion]*(pow(abs(data.AIOMFAC_ION_CHARGE[ion]),2.0)*0.5)*pre_summation6_ion)*c_val*a_val
    for anion, a_val in anion_molality_old.items()
    for cation, c_val in cation_molality_old.items())
    for ion, value in ion_molalities_old.items()
    }
    step=0
    pdb.set_trace()
    summation6_ion_persistent=np.zeros((len(ion_molalities)),)
    for value in ion_molalities:
        summation6_ion_persistent[step]=np.sum(((Cca_persistent*abs(ion_charge[step]))+Cca_dash_persistent*(np.power(abs(ion_charge[step]),2.0)*0.5)*pre_summation6_ion_persistent)*matrix_temp)
        step+=1
    pdb.set_trace()
    #Cation/anion specific calculations

    summation3_ion_cation = {
    cation: sum (Bca[cation][anion]*value for anion, value in anion_molality_old.items())
    for cation, value_c in cation_molality_old.items()
    }
    pdb.set_trace()
    #In the following we need to multiply across rows
    summation3_ion_cation_persistent = np.sum(np.multiply(Bca_persistent,anion_molality[np.newaxis,:]),axis=1)

    summation3_ion_anion = {
    anion: sum (Bca[cation][anion]*value for cation, value in cation_molality_old.items())
    for anion, value_a in anion_molality_old.items()
    }
    pdb.set_trace()
    summation3_ion_anion_persistent = np.sum(np.multiply(Bca_persistent,cation_molality[:,np.newaxis]),axis=0)

    summation3_ion=summation3_ion_cation.copy()
    summation3_ion.update(summation3_ion_anion)
    summation3_ion_persistent =summation3_ion_cation_persistent.copy()
    summation3_ion_persistent=np.append(summation3_ion_persistent,summation3_ion_anion_persistent)
	
    summation5_ion_cation = {
        cation : sum(Cca[cation][anion]*value*pre_summation5_ion for anion, value in anion_molality_old.items())
    for cation, value_c in cation_molality_old.items()
    }
    pdb.set_trace()
    summation5_ion_cation_persistent=np.sum(np.multiply(Cca_persistent,anion_molality[np.newaxis,:])*pre_summation5_ion_persistent,axis=1)
    
    summation5_ion_anion = {
    anion : sum(Cca[cation][anion]*value*pre_summation5_ion for cation, value in cation_molality_old.items())
    for anion, value_a in anion_molality_old.items()
    }
    pdb.set_trace()
    summation5_ion_anion_persistent=np.sum(np.multiply(Cca_persistent,cation_molality[:,np.newaxis])*pre_summation5_ion_persistent,axis=0)

    summation5_ion=summation5_ion_cation.copy()
    summation5_ion.update(summation5_ion_anion)
    summation5_ion_persistent =summation5_ion_cation_persistent.copy()
    summation5_ion_persistent=np.append(summation5_ion_persistent,summation5_ion_anion_persistent)
    pdb.set_trace()

    summation7_ion_cation_1 = {
        cation: sum(Rc_cdash[cation][ion]*value for ion, value in cation_molality_old.items())
    for cation, value_c in cation_molality_old.items()
    }
    pdb.set_trace()
    summation7_ion_cation_1_persistent =np.sum(np.multiply(Rc_cdash_persistent,cation_molality[np.newaxis,:]),axis=1)

    summation7_ion_cation_2 = {
    cation: sum(Qc_cdash[cation][c][a]*c_val*a_val
    for c,c_val in cation_molality_old.items()
    for a,a_val in anion_molality_old.items())
    for cation, value_c in cation_molality_old.items()
    }
    pdb.set_trace()
    summation7_ion_cation_2_persistent=np.zeros((len(cation_molality)),)
    step=0
    for value in cation_molality:
        summation7_ion_cation_2_persistent[step]=np.sum(np.multiply(Qc_cdash_persistent[step,:,:],matrix_temp))
    pdb.set_trace()

    summation7_ion_cation = {
    cation: summation7_ion_cation_1[cation]+summation7_ion_cation_2[cation]
    for cation, value_c in cation_molality_old.items()
    }
    pdb.set_trace()
    summation7_ion_cation_persistent=summation7_ion_cation_1_persistent+summation7_ion_cation_2_persistent

    summation7_ion_anion = {
    anion: sum(Qc_cdash[c][cdash][anion]*c_val*cdash_val
    for c, c_val in cation_molality_old.items()
    for cdash, cdash_val in cation_molality_old.items())
    for anion , a_val in anion_molality_old.items()
    }
    pdb.set_trace()
    summation7_ion_anion_persistent=np.zeros((len(anion_molality)),)
    step=0
    for value in anion_molality:
        summation7_ion_anion_persistent[step]=np.sum(np.multiply(Qc_cdash_persistent[:,:,step],matrix_temp))

    summation7_ion=summation7_ion_cation.copy()
    summation7_ion.update(summation7_ion_anion)
    summation7_ion_persistent=summation7_ion_cation_persistent.copy()
    summation7_ion_persistent=np.append(summation7_ion_persistent,summation7_ion_anion_persistent)

    Ln_gamma_i_MR={
    ion: summation1_ion[ion]+summation2_ion[ion]+summation3_ion[ion]+summation4_ion[ion]+\
        summation5_ion[ion]+summation6_ion[ion]+summation7_ion[ion]
    for ion, value in ion_molalities_old.items()
        if value > 0
    }
    pdb.set_trace()
    Ln_gamma_i_MR_persistent=summation1_ion_persistent+summation2_ion_persistent[:,0]+summation3_ion_persistent+summation4_ion_persistent[:,0]+\
    summation5_ion_persistent+summation6_ion_persistent+summation7_ion_persistent
    pdb.set_trace()
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
    Ln_gamma_s_LR_persistent=((2.0*A*(molw_array_org*1e-3))/(np.power(b,3.0)))*summation_org

    #--Ionic activity coefficients--
    Ln_gamma_i_LR={
        ion: (-1.0*pow((abs(data.AIOMFAC_ION_CHARGE[ion])),2.0)*summation_ion_1)/summation_ion_2 #on the mole fraction scale
        for ion, value in ion_molalities.items()
        }
    Ln_gamma_i_LR_persistent=(-1.0*np.power(abs(ion_charge),2.0)*summation_ion_1)/summation_ion_2

    #Now calculate the conversion factor to change the reference state from mole fraction to molality
    total_abundance_org = sum(value for key, value in organic_compounds.iteritems())
    x_i_org = {
        compound: value/total_abundance_org for compound, value in organic_compounds.items()
        }
    
    x_i_org_persistent=np.array(abundance_array[0:num_orgs,0])/(np.sum(np.array(abundance_array[0:num_orgs,0])))
	
    mean_mw_solvent=sum(x_i_org[compound]*(compound.molwt*1e-3) for compound, value in organic_compounds.items())
    mean_mw_solvent_persistent=np.sum(np.multiply(x_i_org_persistent,molw_array_org)*1e-3)
	
    sum_molalities=sum(value for ion, value in ion_molalities_old.items())
    sum_molalities_persistent=np.sum(ion_molalities)

    ionic_conversion_factor=log((0.01801528/mean_mw_solvent)+0.01801528*sum_molalities)
    ionic_conversion_factor_persistent=log((0.01801528/mean_mw_solvent_persistent)+0.01801528*sum_molalities_persistent)


    Ln_gamma_i_MR_LR= {
        ion: Ln_gamma_i_MR[ion]+Ln_gamma_i_LR[ion]-ionic_conversion_factor
        for ion, value in ion_molalities.items()
        if value > 0
        }
    #Now map this ontp pybel objects as keys rather than integers
    Ln_gamma_i_MR_LR_persistent=Ln_gamma_i_MR_persistent+Ln_gamma_i_LR_persistent-ionic_conversion_factor

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
    pdb.set_trace()
    Ln_gamma_s_MR_LR_persistent = Ln_gamma_s_LR_persistent+Ln_gamma_s_MR_persistent
    Ln_gamma_tot_MR_LR_persistent = np.append(Ln_gamma_s_MR_LR_persistent,Ln_gamma_i_MR_LR_persistent)
        
    # Create a dictionary that holds both ionic and organic values
    
    Ln_gamma_tot_MR_LR={}
    for compound, abundance in organic_compounds.items():
        Ln_gamma_tot_MR_LR[compound] = Ln_gamma_s_MR_LR[compound]
    
    for ion, matches in m_inorg.items():
        Ln_gamma_tot_MR_LR[ion]= Ln_gamma_i_MR_LR_keys[ion]


    #NOTE we only return the organic activity coefficients here. For including inorganic
    #ions, this will also need to be passed
    #pdb.set_trace()
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
        x=c.compound
        c.update(exp(Activity_coefficients_sr[x]+Activity_coefficients_mr_lr[x]))
        step+=1

    #for compound, abundance in m.items(0
    #pdb.set_trace()
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

    #pdb.set_trace()
    #Activity_coefficients_mr_lr=aiomfac_mr(organic_compounds, inorganic_ions, temperature)
    step=0
    for c in m:
        x=c.compound
        c.update(exp(Activity_coefficients_sr[x]))
        
        step+=1

    #for compound, abundance in m.items(0

    return Activity_coefficients_sr, m



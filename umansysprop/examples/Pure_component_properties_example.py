##########################################################################################
#											 #
#    Example file that loads in SMILES strings and then calculates pure component        #
#    properties for a given temperature                                                  # 
#                                                                                        #
#                                                                                        #
#    Copyright (C) 2016  David Topping : david.topping@manchester.ac.uk                  #
#                                      : davetopp80@gmail.com     			 #
#    Personal website: davetoppingsci.com                                                #
#											 #
#    This program is free software: you can redistribute it and/or modify                #
#    it under the terms of the GNU Affero General Public License as published            #
#    by the Free Software Foundation, either version 3 of the License, or                #
#    (at your option) any later version.                                                 # 
#                                                                                        #   
#    This program is distributed in the hope that it will be useful,                     #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                      # 
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                       #
#    GNU Affero General Public License for more details.                                 #
#                                                                                        #
#    You should have received a copy of the GNU Affero General Public License            # 
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.               #
#                                                                                        #
#                                                                                        #
##########################################################################################

# This file demonstrates how to call existing routines in UManSysProp to extract pure
# component properties. Please note this file is not designed to be run directly from
# this sub-directory of the main distribution. It is deisgned to illustrate how one might
# call the seperate methods according to your UManSysProp location.
#
# The conversion from uploaded SMILES to a Pybel object is also different from the mechanism
# used in the files in the 'tools' directory. Those files are designed to work with
# forms on the current website. The principle is the same, we simply call Pybel directly.
#
# Current properties calculated:
# - Pure component vapour pressure
# - Boiling point
# - Subcooled liquid density
# - Molecular weight and chemical formula
#
# Last modification 21/12/16

import glob
from openbabel import pybel
import collections
import sys
# Here you will need to put the relevant path in for your UManSysProp distribution. Mine is
# given as an example  - change this!
sys.path.append('/Users/user/Documents/GitHub/UManSysProp_public/')
from umansysprop import boiling_points
from umansysprop import vapour_pressures
from umansysprop import critical_properties
from umansysprop import liquid_densities
import pdb

##########################################################################################
# 1. Read in the property data, seperating only .prop files from the rest.
# I have defined the .prop files myself, just for ease of use for this project.
# It is just a textfile with compound name/reference as first column, SMILES as second.

onlyfiles = [f for f in glob.glob('*.prop')]

#1) extract the data from each file and start generating a dictionary of Pybel objects

step=0
filenames=[]
Compound_reference=[]
smiles_array=[]
property_array=[]
Pybel_object_dict=dict()

for filename in onlyfiles: # If you have more than one file, for whatever reason

   SMILES_flag=0
   filenames.append(filename[:])
   text=open(filename[:],'r')
   
   for line in text.readlines():
       input = line.split()
       # Keep a list of the information
       Compound_reference.append(input[0])
       smiles_array.append(input[1])
       # Now create Pybel objects which are used in all property predictive techniquexs
       Pybel_object=pybel.readstring('smi',input[1])
       Pybel_object_dict[input[1]]=Pybel_object
       
##########################################################################################       
# 2) Create a dictionary of properties based on these Pybel objects

# NOTE: For some of the vapour pressure values, you need to perform a boiling point estimation first
# It is therefore wise to do this initially

# 2a) Boiling points [(K)]

boiling_point_dict=collections.defaultdict(lambda: collections.defaultdict())

for smiles in smiles_array:

    boiling_point_dict[smiles]['joback_and_reid']=boiling_points.joback_and_reid(Pybel_object_dict[smiles])
    boiling_point_dict[smiles]['stein_and_brown']=boiling_points.stein_and_brown(Pybel_object_dict[smiles])
    boiling_point_dict[smiles]['nannoolal']=boiling_points.nannoolal(Pybel_object_dict[smiles])
    

# 2b) Vapour pressures [log10 (atm) at a specific temperature]
# For those vapour pressure methods that require a boiling point, we have 3D dictionaries

vapour_pressure_dict_BP=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
vapour_pressure_dict=collections.defaultdict(lambda: collections.defaultdict())
temperature=298.15

for smiles in smiles_array:

    vapour_pressure_dict_BP[smiles]['VP_Nannoolal']['BP_Nannoolal']=vapour_pressures.nannoolal(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['nannoolal'])
    vapour_pressure_dict_BP[smiles]['VP_Nannoolal']['BP_Stein_Brown']=vapour_pressures.nannoolal(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['stein_and_brown'])
    vapour_pressure_dict_BP[smiles]['VP_Nannoolal']['BP_Joback_Reid']=vapour_pressures.nannoolal(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['joback_and_reid'])
    
    vapour_pressure_dict_BP[smiles]['VP_Myrdal_Yalkowsky']['BP_Nannoolal']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['nannoolal'])
    vapour_pressure_dict_BP[smiles]['VP_Myrdal_Yalkowsky']['BP_Stein_Brown']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['stein_and_brown'])
    vapour_pressure_dict_BP[smiles]['VP_Myrdal_Yalkowsky']['BP_Joback_Reid']=vapour_pressures.myrdal_and_yalkowsky(Pybel_object_dict[smiles], temperature, boiling_point_dict[smiles]['joback_and_reid'])
    
    vapour_pressure_dict[smiles]['EVAP']=vapour_pressures.evaporation(Pybel_object_dict[smiles], temperature)


# 2c) Subooled liquid density [(g/cc) at a specific temperature]
# Some require a pre-calculated critical properties and thus boiling points which are fed into the function directly

density_CP=collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
density=collections.defaultdict(lambda: collections.defaultdict())

for smiles in smiles_array:

    compound=Pybel_object_dict[smiles] #Using this simply to try and keep line extensions short as possible
    boiling_point=boiling_point_dict[smiles]['nannoolal'] #Selecting this as an example

    density[smiles]['Girolami']=liquid_densities.girolami(compound)
    
    density_CP[smiles]['Schroeder']['CP_Nannoolal']=liquid_densities.schroeder(compound, temperature, critical_properties.nannoolal(compound, boiling_point))
    density_CP[smiles]['Schroeder']['CP_Joback_Reid']=liquid_densities.schroeder(compound, temperature, critical_properties.nannoolal(compound, boiling_point))
    density_CP[smiles]['Le_Bas']['CP_Nannoolal']=liquid_densities.le_bas(compound, temperature, critical_properties.nannoolal(compound, boiling_point))
    density_CP[smiles]['Le_Bas']['CP_Joback_Reid']=liquid_densities.le_bas(compound, temperature, critical_properties.nannoolal(compound, boiling_point))
    density_CP[smiles]['Tyn_Calus']['CP_Nannoolal']=liquid_densities.tyn_and_calus(compound, temperature, critical_properties.nannoolal(compound, boiling_point))
    density_CP[smiles]['Tyn_Calus']['CP_Joback_Reid']=liquid_densities.tyn_and_calus(compound, temperature, critical_properties.nannoolal(compound, boiling_point))

# 2d) Pybel object properties. By creating a pybel object you can extract a range of markers/features of the molecule
# In Python, these features of the Pybel object class can be found by typing dir(inset name of your pybel object here)
# For example:

general_dict=collections.defaultdict(lambda: collections.defaultdict())

for smiles in smiles_array:

    general_dict[smiles]['Mw']=Pybel_object_dict[smiles].molwt
    general_dict[smiles]['Formula']=Pybel_object_dict[smiles].formula
    
##########################################################################################
# 3) Example output
# The debugging stop point below allows you to interogate some derived values.
# For example, in the command line type:
# density_CP[smiles_array[4]]['Le_Bas']['CP_Nannoolal']
# vapour_pressure_dict_BP[smiles_array[4]]['VP_Myrdal_Yalkowsky']['BP_Nannoolal']
# To extract the density of the 5th compound in the smiles_array according to:
# - The method of Le_Bas using critical properties by Nannoolal
# and to extract the vapour pressure of the 5th compound in the smiles_array according to:
# - The method of Myrdal and Yalkowsky, using the boiling point estimation by Nannoolal
#
# To finish the program, simply type 'c' to continue. The following code then saves some
# of the vapour pressure estimations to a file. 
    
pdb.set_trace()

file_name = open('Example_output.txt','w')
file_name.write(str('SMILES')+'\t'+str('VP_Nannoolal_BP_Nannoolal')+'\t'+str('EVAP')) 
    
for smiles in smiles_array:

    file_name.write('\n')
    file_name.write(str(smiles)+'\t'+str(vapour_pressure_dict_BP[smiles]['VP_Nannoolal']['BP_Nannoolal'])+'\t'+str(vapour_pressure_dict[smiles]['EVAP']))  

file_name.close()




#    Copyright 2015 Adriana Supady 
#    adriana.supady@gmail.com
#
#    This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.
''' This module handles the molecule and its 3D structures '''

import pybel
import openbabel
import numpy as np
import glob
import os
import shutil
import subprocess
from copy import deepcopy
from random import randint, choice

from new_molecule_parametrize import parametrize, template_sdf
from angle import angle_measure, angle_set
from check_geometry import is_geometry_valid
from genetic_operations import crossover_vec, mutation_tor_vec, mutation_cistrans_vec
from pyaims import AimsObject
from pyff import FFObject
from utilities import *



class MoleculeDescription:
    newline = "NEWLINE"
    def __init__(self, parameter_file=None, **kwargs):
        if parameter_file != None:
            with open(parameter_file) as fin:
                parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
            self.smile = parameter_dict['smile']
            self.custom = boolean(parameter_dict['custom'])
            self.rmsd_cutoff_uniq = float(parameter_dict['rmsd_cutoff_uniq'])
            self.distance_cutoff_1 = float(parameter_dict['distance_cutoff_1'])
            self.distance_cutoff_2 = float(parameter_dict['distance_cutoff_2'])
            self.chiral = boolean(parameter_dict['chiral'])
            self.rmsd_type = parameter_dict['rmsd_type']
            self.smart_torsion = pybel.Smarts(parameter_dict['smart_torsion']) 
            self.smart_cistrans = pybel.Smarts(parameter_dict['smart_cistrans'])
            self.smart_custom = pybel.Smarts(parameter_dict['smart_custom'])
        else:
            for key in kwargs.keys(): 
                if key not in ["template_sdf_string"]:
                    setattr(self, str(key), kwargs[key])
                else:
                    setattr(self, str(key), kwargs[key].replace(MoleculeDescription.newline, "\n")) 
                
            
    def __repr__(self):
        repr_list = []
        for att_name in self.__dict__.keys():
            if att_name in ["smile",  "rmsd_type", "smart_torsion", "smart_cistrans", "smart_custom"]:
                repr_list.append('%s="%s"' %(att_name,getattr(self, att_name)))
            elif att_name in ["rmsd_cutoff_uniq","distance_cutoff_1", "distance_cutoff_2", "atoms", "bonds", "custom", "chiral", "torsion", "cistrans", "custom_torsion"]:
                repr_list.append('%s=%s' %(att_name,repr(getattr(self, att_name))))
            elif att_name in ["template_sdf_string"]:
                repr_list.append("%s='%s'" %(att_name,getattr(self, att_name).replace("\n", MoleculeDescription.newline,)))
            elif att_name in ["matrix_of_connections"]:
                repr_list.append('%s=%s' %(att_name,' '.join(repr(getattr(self, att_name)).split())))
            else: 
                "Unknown type of attribute"
        return "%s(%s)" %(self.__class__.__name__,', '.join(repr_list) )
   

    def get_mol_parameters(self):
        (self.atoms, self.bonds, self.torsion, self.cistrans,  self.custom_torsion) = parametrize(self.smile, self.custom, self.smart_torsion, self.smart_cistrans, self.smart_custom)
        if self.custom:
            self.torsion = self.custom_torsion
    def create_template_sdf(self):
        (self.matrix_of_connections, self.template_sdf_string)=template_sdf(self.smile, self.atoms, self.bonds, self.distance_cutoff_1, self.distance_cutoff_2)  
        
class Structure:
    index = 0
    newline = "NEWLINE"
    def __init__(self, arg=None, **kwargs):
        if isinstance(arg, MoleculeDescription):
            self.mol_info = arg
            Structure.index += 1
            self.index = Structure.index
            
        elif isinstance(arg,Structure):
            self.mol_info = arg.mol_info 
            Structure.index += 1
            self.index = Structure.index
            for att_name in arg.__dict__.keys():
                if att_name != "mol_info" and att_name != "index":
                    setattr(self, str(att_name),deepcopy(getattr(arg, str(att_name))))
        else:            
            "The initialization cannot be performed. Check your input"
        for key in kwargs.keys(): 
            if key not in ["xyz_string", "initial_xyz_string"]:
                setattr(self, str(key), kwargs[key])
            else:
                setattr(self, str(key), kwargs[key].replace(Structure.newline, "\n")) 
                    
                
    def __repr__(self):
         

        repr_list = []
        for att_name in self.__dict__.keys():
            
            if att_name in ["energy","index", "values_cistrans", "values_tor", "initial_values_cistrans", "initial_values_tor"]:
                repr_list.append('%s=%s' %(att_name,repr(getattr(self, att_name))))
            elif att_name in ["xyz_string", "initial_xyz_string"]:
                repr_list.append("%s='%s'" %(att_name,getattr(self, att_name).replace("\n", Structure.newline)))
     
            else: 
                "Unknown type of attribute"
        return "%s(mol, %s)" %(self.__class__.__name__,', '.join(repr_list) )

    def __str__(self):
        return "%s %d"%(self.__class__.__name__,self.index)
      
    def __float__(self):
        return float(self.energy)
        
    def generate_random_structure(self, cistrans1_value, cistrans2_value):
        
        self.values_cistrans, self.values_tor = [],[]
    
        for x in range(len(self.mol_info.cistrans)):
            self.values_cistrans.append(float(choice([cistrans1_value,cistrans2_value])))
        for y in range(len(self.mol_info.torsion)):
            self.values_tor.append(randint(0.0, 359.0)-179.0) 
        
        self.xyz_string=angle_set(self.mol_info.template_sdf_string, self.mol_info.cistrans, self.mol_info.torsion, self.values_cistrans, self.values_tor)
        self.values_cistrans, self.values_tor=angle_measure(self.xyz_string, self.mol_info.cistrans, self.mol_info.torsion)
        
    def generate_structure_from_values(self, values_cis, values_tor):
        
        self.values_cistrans, self.values_tor = values_cis, values_tor
        self.xyz_string=angle_set(self.mol_info.template_sdf_string, self.mol_info.cistrans, self.mol_info.torsion,self.values_cistrans, self.values_tor)
        self.values_cistrans, self.values_tor=angle_measure(self.xyz_string, self.mol_info.cistrans, self.mol_info.torsion)
        
    def is_geometry_valid(self):
        return is_geometry_valid(self.xyz_string, self.mol_info.atoms, self.mol_info.bonds, self.mol_info.matrix_of_connections, self.mol_info.distance_cutoff_1, self.mol_info.distance_cutoff_2)

    def __eq__(self, other):
        if self.mol_info.rmsd_type == 'cartesian':
            obj1, obj2 = self, other
            
            if hasattr(self,"initial_xyz_string"):
                obj1, obj2 = obj2, obj1           
            if hasattr(obj1,"initial_xyz_string"):
                raise Exception("Both structures are already relaxed, no need for checking the uniquness)")
            smart = obj1.mol_info.smile

            filename1="ref.xyz"
            with open(filename1, 'w') as f:
                f.write(obj1.xyz_string)
            filename2="target.xyz"
            with open(filename2, 'w') as f:
                f.write(obj2.xyz_string)
            a=subprocess.Popen("obfit "+str("\"")+str(smart)+str("\"")+str(" ")+str(filename1)+str(" ")+str(filename2), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (output, err)=a.communicate()
            try:
                rmsd1=float(err.split('\n')[0].split(' ')[1])
            except:
                pass
            rmsd = rmsd1
            
            if hasattr(obj2,"initial_xyz_string"):
                filename3="target_ini.xyz"
                with open(filename3, 'w') as f:
                    f.write(obj2.initial_xyz_string)
                a=subprocess.Popen("obfit "+str("\"")+str(smart)+str("\"")+str(" ")+str(filename1)+str(" ")+str(filename3), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (output, err)=a.communicate()
                try:
                    rmsd2=float(err.split('\n')[0].split(' ')[1])
                    if rmsd2<rmsd1: rmsd=rmsd2
                except:
                    pass
                
            
            
            if not obj1.mol_info.chiral:
                filename2mirr = "target_mirror.xyz"
                mirrorer(filename2, filename2mirr) 
                a=subprocess.Popen("obfit "+str("\"")+str(smart)+str("\"")+str(" ")+str(filename1)+str(" ")+str(filename2mirr), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (output, err)=a.communicate()
                try:
                    rmsd_mirr1=float(err.split('\n')[0].split(' ')[1])
                except:
                    pass
                rmsd_mirr=rmsd_mirr1
                if hasattr(obj2, "initial_xyz_string"):
                    filename3mirr = "target_ini_mirror.xyz" 
                    mirrorer(filename3, filename3mirr)
                    a=subprocess.Popen("obfit "+str("\"")+str(smart)+str("\"")+str(" ")+str(filename1)+str(" ")+str(filename3mirr), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    (output, err)=a.communicate()
                    try:
                        rmsd_mirr2=float(err.split('\n')[0].split(' ')[1])
                        if rmsd_mirr2<rmsd_mirr1: rmsd_mirr = rmsd_mirr2
                    except:
                        pass
                    remover_file(filename3mirr)
                if rmsd_mirr<rmsd: rmsd = rmsd_mirr
                
                remover_file(filename2mirr)
                
                
            for file in filename1, filename2:
                remover_file(file)

            remover_file(filename3)
            if rmsd > obj1.mol_info.rmsd_cutoff_uniq: return False
            else: return True     
           
   
        if self.mol_info.rmsd_type == 'torsional':
            obj1, obj2 = self, other
            
            if hasattr(self,"initial_xyz_string"):
                obj1, obj2 = obj2, obj1           
            if hasattr(obj1,"initial_xyz_string"):
                raise Exception("Both structures are already relaxed, no need for checking the uniquness)")
            
            
            rmsd1=tor_rmsd(2,get_vec(np.concatenate((obj1.values_cistrans,obj1.values_tor), axis=0), np.concatenate((obj2.values_cistrans,obj2.values_tor), axis=0)))
            rmsd = rmsd1
            if hasattr(obj2,"initial_xyz_string"):
                rmsd2=tor_rmsd(2,get_vec(np.concatenate((obj1.values_cistrans,obj1.values_tor), axis=0), np.concatenate((obj2.initial_values_cistrans,obj2.initial_values_tor), axis=0)))
                if rmsd2<rmsd1: rmsd=rmsd2
  
             
            if not obj1.mol_info.chiral:
                rmsd_mirr1=tor_rmsd(2,get_vec(np.concatenate((obj1.values_cistrans,obj1.values_tor), axis=0), -1*np.concatenate((obj2.values_cistrans,obj2.values_tor), axis=0))) 
                rmsd_mirr2=tor_rmsd(2,get_vec(np.concatenate((obj1.values_cistrans,obj1.values_tor), axis=0), -1*np.concatenate((obj2.initial_values_cistrans,obj2.initial_values_tor), axis=0)))
                rmsd_mirr = min(rmsd_mirr1, rmsd_mirr2)
                if rmsd_mirr < rmsd : rmsd=rmsd_mirr
            
            
            if rmsd > obj1.mol_info.rmsd_cutoff_uniq: return False
            else: return True
        
    def send_to_blacklist(self,address,array):
        try: 
            os.mkdir(address)
        except OSError:
            pass
        array.append(self)
        cnt_black = len(glob.glob(address+"/*.xyz"))
        filename="tmp.xyz"
        with open(filename, 'w') as f:
            f.write(self.xyz_string)
        shutil.copy("tmp.xyz", str(address)+"/black_"+str(cnt_black)+".xyz")

        if not self.mol_info.chiral:
            new_file_name = "mirror.xyz" 
            mirrorer("tmp.xyz", new_file_name)
            cnt_black = len(glob.glob(address+"/*.xyz"))
            shutil.copy("mirror.xyz", str(address)+"/black_"+str(cnt_black)+".xyz")
            os.remove("mirror.xyz")
        remover_file("tmp.xyz")
 
        
    def perform_aims(self, parameter_file, execution_string, dirname):
        aims_object = AimsObject(parameter_file)
        aims_object.generate_input(self.xyz_string)
        aims_object.build_storage(dirname)
        aims_object.run_aims(execution_string)
        aims_object.clean_and_store()
        self.energy = aims_object.get_energy()
         
        self.initial_xyz_string = self.xyz_string
        self.initial_values_cistrans, self.initial_values_tor  = self.values_cistrans, self.values_tor

        self.xyz_string=fhiaims_to_xyz_string(aims_object.get_aims_string_opt())
        self.values_cistrans, self.values_tor=angle_measure(self.xyz_string, self.mol_info.cistrans, self.mol_info.torsion)
        
    def perform_ff(self, parameter_file):
        ff_object  = FFObject(parameter_file)
        ff_object.run_ff(self.xyz_string)
        self.energy = ff_object.get_energy()
        
        self.initial_xyz_string = self.xyz_string
        self.initial_values_cistrans, self.initial_values_tor  = self.values_cistrans, self.values_tor
        
        self.xyz_string = ff_object.get_xyz_string_opt()
        self.values_cistrans, self.values_tor=angle_measure(self.xyz_string, self.mol_info.cistrans, self.mol_info.torsion)
        
    def __cmp__(self, other):
        return cmp(self.energy, other.energy)
    
    def crossover(self, other):
        start1, start2 = crossover_vec(self.values_cistrans, other.values_cistrans)
        end1, end2 = crossover_vec(self.values_tor, other.values_tor)
              
        child1_vec=np.append(start1,end1)
        child2_vec=np.append(start2,end2)
        
        child1=Structure(self.mol_info)
        child1.generate_structure_from_values(start1, end1)
        child2=Structure(self.mol_info)
        child2.generate_structure_from_values(start2, end2)
        
        return (child1, child2)

    
    def mutation_tor(self,max_mutations_torsions):
        dt = len(self.mol_info.torsion)
        if dt > 0:
            mut_vec=mutation_tor_vec(self.values_tor,max_mutations_torsions)
            self.generate_structure_from_values(self.values_cistrans, mut_vec )

    def mutation_cistrans(self,max_mutations_cistrans):
        dc = len(self.mol_info.cistrans)
        if dc > 0:
            mut_vec=mutation_cistrans_vec(self.values_cistrans,max_mutations_cistrans)
            self.generate_structure_from_values(mut_vec, self.values_tor)

            
                
    
from __future__ import print_function
from __future__ import absolute_import
import shutil
import uuid
import getpass
import socket
import os
import sys
import numpy as np
import scipy.spatial as scsp
import copy
import time
import subprocess
import shlex
import yaml
from scipy.spatial.distance import cdist

from morfeus import BuriedVolume
from morfeus import Pyramidalization
from morfeus import ConeAngle
from morfeus import SolidAngle
from morfeus import Sterimol
from morfeus import SASA
from morfeus import Dispersion

import scipy.spatial as scsp
import scipy.linalg as scli
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi



kcal_to_eV=0.0433641153
kB=8.6173303e-5 #eV/K
T=298.15
kBT=kB*T
AToBohr=1.889725989

# Generates the adjacency matrix based on UFF bond radii
# Inputs:       Elements: N-element List strings for each atom type
#               Geometry: Nx3 np.array holding the geometry of the molecule
#               File:  Optional. If Table_generator encounters a problem then it is often useful to have the name of the file the geometry came from printed. 
def Table_generator(Elements,Geometry,File=None,Radii_dict=False):
    module_path =  '/'.join(os.path.abspath(__file__).split('/')[:-1])
    # Initialize UFF bond radii (Rappe et al. JACS 1992)
    # NOTE: Units of angstroms 
    # NOTE: These radii neglect the bond-order and electronegativity corrections in the original paper. Where several values exist for the same atom, the largest was used. 
    Radii = {  'H':0.354, 'He':0.849,\
              'Li':1.336, 'Be':1.074,                                                                                                                          'B':0.838,  'C':0.757,  'N':0.700,  'O':0.658,  'F':0.668, 'Ne':0.920,\
              'Na':1.539, 'Mg':1.421,                                                                                                                         'Al':1.244, 'Si':1.117,  'P':1.117,  'S':1.064, 'Cl':1.044, 'Ar':1.032,\
               'K':1.953, 'Ca':1.761, 'Sc':1.513, 'Ti':1.412,  'V':1.402, 'Cr':1.345, 'Mn':1.382, 'Fe':1.335, 'Co':1.241, 'Ni':1.164, 'Cu':1.302, 'Zn':1.193, 'Ga':1.260, 'Ge':1.197, 'As':1.211, 'Se':1.190, 'Br':1.192, 'Kr':1.147,\
              'Rb':2.260, 'Sr':2.052,  'Y':1.698, 'Zr':1.564, 'Nb':1.473, 'Mo':1.484, 'Tc':1.322, 'Ru':1.478, 'Rh':1.332, 'Pd':1.338, 'Ag':1.386, 'Cd':1.403, 'In':1.459, 'Sn':1.398, 'Sb':1.407, 'Te':1.386,  'I':1.382, 'Xe':1.267,\
              'Cs':2.570, 'Ba':2.277, 'La':1.943, 'Hf':1.611, 'Ta':1.511,  'W':1.526, 'Re':1.372, 'Os':1.372, 'Ir':1.371, 'Pt':1.364, 'Au':1.262, 'Hg':1.340, 'Tl':1.518, 'Pb':1.459, 'Bi':1.512, 'Po':1.500, 'At':1.545, 'Rn':1.42,\
              'default' : 0.7 }

    # SAME AS ABOVE BUT WITH A SMALLER VALUE FOR THE Al RADIUS ( I think that it tends to predict a bond where none are expected
    Radii = {  'H':0.39, 'He':0.849,\
              'Li':1.336, 'Be':1.074,                                                                                                                          'B':0.838,  'C':0.757,  'N':0.700,  'O':0.658,  'F':0.668, 'Ne':0.920,\
              'Na':1.539, 'Mg':1.421,                                                                                                                         'Al':1.15,  'Si':1.050,  'P':1.117,  'S':1.064, 'Cl':1.044, 'Ar':1.032,\
               'K':1.953, 'Ca':1.761, 'Sc':1.513, 'Ti':1.412,  'V':1.402, 'Cr':1.345, 'Mn':1.382, 'Fe':1.335, 'Co':1.241, 'Ni':1.164, 'Cu':1.302, 'Zn':1.193, 'Ga':1.260, 'Ge':1.197, 'As':1.211, 'Se':1.190, 'Br':1.192, 'Kr':1.147,\
              'Rb':2.260, 'Sr':2.052,  'Y':1.698, 'Zr':1.564, 'Nb':1.473, 'Mo':1.484, 'Tc':1.322, 'Ru':1.478, 'Rh':1.332, 'Pd':1.338, 'Ag':1.386, 'Cd':1.403, 'In':1.459, 'Sn':1.398, 'Sb':1.407, 'Te':1.386,  'I':1.382, 'Xe':1.267,\
              'Cs':2.570, 'Ba':2.277, 'La':1.943, 'Hf':1.611, 'Ta':1.511,  'W':1.526, 'Re':1.372, 'Os':1.372, 'Ir':1.371, 'Pt':1.364, 'Au':1.262, 'Hg':1.340, 'Tl':1.518, 'Pb':1.459, 'Bi':1.512, 'Po':1.500, 'At':1.545, 'Rn':1.42,\
              'default' : 0.7 }

    # Use Radii json file in Lib folder if sepcified
    if Radii_dict == True:
      if os.path.isfile(module_path+'/Radii.json') is False:
         print("ERROR: {}/Radii.json doesn't exist, check for Radii.json file in the Library".format(module_path))
         quit()
      Radii = read_alljson(module_path+'/Radii.json')


    Max_Bonds = {  'H':2,    'He':1,\
                  'Li':None, 'Be':None,                                                                                                                'B':4,     'C':4,     'N':4,     'O':2,     'F':1,    'Ne':1,\
                  'Na':None, 'Mg':None,                                                                                                               'Al':4,    'Si':4,  'P':None,  'S':None, 'Cl':1,    'Ar':1,\
                   'K':None, 'Ca':None, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':3,    'Ge':None, 'As':None, 'Se':None, 'Br':1,    'Kr':None,\
                  'Rb':None, 'Sr':None,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':1,    'Xe':None,\
                  'Cs':None, 'Ba':None, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }
                     
    # Scale factor is used for determining the bonding threshold. 1.2 is a heuristic that give some lattitude in defining bonds since the UFF radii correspond to equilibrium lengths. 
    scale_factor = 1.2

    # Print warning for uncoded elements.
    for i in Elements:
        if i not in Radii.keys():
            print( "ERROR in Table_generator: The geometry contains an element ({}) that the Table_generator function doesn't have bonding information for. This needs to be directly added to the Radii".format(i)+\
                  " dictionary before proceeding. Exiting...")
            quit()

    # Generate distance matrix holding atom-atom separations (only save upper right)
    Dist_Mat = np.triu(cdist(Geometry,Geometry))
    
    # Find plausible connections
    x_ind,y_ind = np.where( (Dist_Mat > 0.0) & (Dist_Mat < max([ Radii[i]**2.0 for i in Radii.keys() ])) )

    # Initialize Adjacency Matrix
    Adj_mat = np.zeros([len(Geometry),len(Geometry)])

    # Iterate over plausible connections and determine actual connections
    for count,i in enumerate(x_ind):
        
        # Assign connection if the ij separation is less than the UFF-sigma value times the scaling factor
        if Dist_Mat[i,y_ind[count]] < (Radii[Elements[i]]+Radii[Elements[y_ind[count]]])*scale_factor:            
            Adj_mat[i,y_ind[count]]=1
    
        if Elements[i] == 'H' and Elements[y_ind[count]] == 'H':
            if Dist_Mat[i,y_ind[count]] < (Radii[Elements[i]]+Radii[Elements[y_ind[count]]])*1.5:
                Adj_mat[i,y_ind[count]]=1

    # Hermitize Adj_mat
    Adj_mat=Adj_mat + Adj_mat.transpose()

    # Perform some simple checks on bonding to catch errors
    problem_dict = { i:0 for i in Radii.keys() }
    conditions = { "H":1, "C":4, "F":1, "Cl":1, "Br":1, "I":1, "O":2, "N":4, "B":4 }
    for count_i,i in enumerate(Adj_mat):

        if Max_Bonds[Elements[count_i]] is not None and sum(i) > Max_Bonds[Elements[count_i]]:
            problem_dict[Elements[count_i]] += 1
            cons = sorted([ (Dist_Mat[count_i,count_j],count_j) if count_j > count_i else (Dist_Mat[count_j,count_i],count_j) for count_j,j in enumerate(i) if j == 1 ])[::-1]
            while sum(Adj_mat[count_i]) > Max_Bonds[Elements[count_i]]:
                sep,idx = cons.pop(0)
                Adj_mat[count_i,idx] = 0
                Adj_mat[idx,count_i] = 0
#        if Elements[count_i] in conditions.keys():
#            if sum(i) > conditions[Elements[count_i]]:


    # Print warning messages for obviously suspicious bonding motifs.
    if sum( [ problem_dict[i] for i in problem_dict.keys() ] ) > 0:
        print( "Table Generation Warnings:")
        for i in sorted(problem_dict.keys()):
            if problem_dict[i] > 0:
                if File is None:
                    if i == "H": print( "WARNING in Table_generator: {} hydrogen(s) have more than one bond.".format(problem_dict[i]))
                    if i == "C": print( "WARNING in Table_generator: {} carbon(s) have more than four bonds.".format(problem_dict[i]))
                    if i == "Si": print( "WARNING in Table_generator: {} silicons(s) have more than four bonds.".format(problem_dict[i]))
                    if i == "F": print( "WARNING in Table_generator: {} fluorine(s) have more than one bond.".format(problem_dict[i]))
                    if i == "Cl": print( "WARNING in Table_generator: {} chlorine(s) have more than one bond.".format(problem_dict[i]))
                    if i == "Br": print( "WARNING in Table_generator: {} bromine(s) have more than one bond.".format(problem_dict[i]))
                    if i == "I": print( "WARNING in Table_generator: {} iodine(s) have more than one bond.".format(problem_dict[i]))
                    if i == "O": print( "WARNING in Table_generator: {} oxygen(s) have more than two bonds.".format(problem_dict[i]))
                    if i == "N": print( "WARNING in Table_generator: {} nitrogen(s) have more than four bonds.".format(problem_dict[i]))
                    if i == "B": print( "WARNING in Table_generator: {} bromine(s) have more than four bonds.".format(problem_dict[i]))
                else:
                    if i == "H": print( "WARNING in Table_generator: parsing {}, {} hydrogen(s) have more than one bond.".format(File,problem_dict[i]))
                    if i == "C": print( "WARNING in Table_generator: parsing {}, {} carbon(s) have more than four bonds.".format(File,problem_dict[i]))
                    if i == "Si": print( "WARNING in Table_generator: parsing {}, {} silicons(s) have more than four bonds.".format(File,problem_dict[i]))
                    if i == "F": print( "WARNING in Table_generator: parsing {}, {} fluorine(s) have more than one bond.".format(File,problem_dict[i]))
                    if i == "Cl": print( "WARNING in Table_generator: parsing {}, {} chlorine(s) have more than one bond.".format(File,problem_dict[i]))
                    if i == "Br": print( "WARNING in Table_generator: parsing {}, {} bromine(s) have more than one bond.".format(File,problem_dict[i]))
                    if i == "I": print( "WARNING in Table_generator: parsing {}, {} iodine(s) have more than one bond.".format(File,problem_dict[i]))
                    if i == "O": print( "WARNING in Table_generator: parsing {}, {} oxygen(s) have more than two bonds.".format(File,problem_dict[i]))
                    if i == "N": print( "WARNING in Table_generator: parsing {}, {} nitrogen(s) have more than four bonds.".format(File,problem_dict[i]))
                    if i == "B": print( "WARNING in Table_generator: parsing {}, {} bromine(s) have more than four bonds.".format(File,problem_dict[i]))
        print( "")

    return Adj_mat

### graph traversal algorithms --- Depth-First Search (DFS), given start atom idx and adjmat to identify all atoms connect to the seed
def DFS(adjmat, start):
    visited = set()
    stack = [start]

    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            neighbors = [i for i, v in enumerate(adjmat[node]) if v]
            stack.extend(neighbors)
    return visited

def sanitize_smiles(smi):
    return mol2smi(smi2mol(smi, sanitize=True), isomericSmiles=False, canonical=True)

def read_xtb_log1(moldir,mol_len,idx_dict):
    os.chdir(moldir)

    metal_idx = idx_dict['metal_idx']
    ligand_indices = idx_dict['ligand_indices']
    br_idx = idx_dict['br_idx']
    
    if not os.path.exists("xtb.log"):
        return(False, None, None, None, None, None, None, None, None, None, None, None)

    read_mul=False
    read_wil=False
    read_dip=False
    muls=[]
    alphas=[]
    covCN = []
    c6 = []
    muls = []
    dip=[0.0,0.0,0.0]
    HOMO_LUMO_gap=None
    occ_energies=[]
    virt_energies=[]
    occ_done=False
    read_orbital_energies=False
    lines_f = []

    with open("xtb.log","r") as f:
        for lc,line in enumerate(f):
            lines_f += [line]

    # get the lc idx for covCN (q, c6aa and alphas)
    for lc,line in enumerate(lines_f):
        if "Z" in line and "covCN" in line and "q" in line and "C6AA" in line and "α(0)" in line:
            mul_end = lc+mol_len
            print ("mul_end={}".format(mul_end))
            break

    # get the lc idx for wbo
    wbo_start = None
    wbo_end = None
    for lc,line in enumerate(lines_f):
        if 'Z sym  total        # sym  WBO       # sym  WBO       # sym  WBO' in line:
            wbo_start = lc+2
            continue
        elif 'Topologies differ in total number of bonds' in line:
            wbo_end = lc-3
            break

    wbo_data = {}
    wbo_values = {'metal_total': 0, 'metal-ligand': 0, 'metal-br': 0}
    if wbo_start is None or wbo_end is None:
        print("Start or End index not found")
        return wbo_values

    else:
        for lc,line in enumerate(lines_f):
            fields = line.split()
            if lc>= wbo_start and lc<= wbo_end:
                if '--' in line and fields[0].isdigit():
                    key_atom_idx = fields[0]
                    key_element = fields[2]
                    key = key_element+'_'+key_atom_idx
                    elements = fields[4:]
                    wbo_data[key] = {}
                    wbo_data[key]['total'] = fields[3]

                    while elements:
                        if elements and elements[0] == '--':
                            elements.pop(0)
                            continue

                        atom_idx = int(elements.pop(0))
                        element = elements.pop(0)
                        value = float(elements.pop(0))
                        wbo_data[key][str(element)+'_'+str(atom_idx)] = value

                elif '--' not in line and fields[0].isdigit():
                    elements = fields[0:]
                    while elements:
                        atom_idx = int(elements.pop(0))
                        element = elements.pop(0)
                        value = float(elements.pop(0))
                        wbo_data[key][str(element)+'_'+str(atom_idx)] = value

        #print (wbo_data)
    
        # parse the wbo of metal-br, metal-ligand, iterate over the wbo dictionary
        for key, value in wbo_data.items():
            # Check if the key contains the metal
            if metal_idx == int(key.split('_')[-1]) - 1: #0 indices
                for value_i, wbo_value in value.items():
                    if value_i == 'total':
                        wbo_values['metal_total'] = float(wbo_value)
                    elif br_idx == int(value_i.split('_')[-1]) - 1:
                        wbo_values['metal-br'] = wbo_value
                    elif int(value_i.split('_')[-1]) - 1 in ligand_indices:
                        wbo_values['metal-ligand']+=wbo_value
        #print (wbo_values)

    # parse other information
    for lc,line in enumerate(lines_f):
        if "convergence criteria cannot be satisfied within" in line:
            break
        if read_orbital_energies and "HL-Gap" in line:
            read_orbital_energies=False
        if read_dip and "molecular quadrupole" in line:
            read_dip=False
        if lc>mul_end-mol_len and lc<=mul_end:
            #print ('mul_line={}'.format(line))
            muls.append(float(line.split()[4])) #q, charge compare to neutral
            covCN.append(float(line.split()[3])) #covCN, coordination number, double bond also count as 1 (different with total wbo: double bound count as 2)
            c6.append(float(line.split()[5])) #C6 dispersion coefficient
            alphas.append(float(line.split()[6])) #α(0)
        if read_dip and "full:" in line and len(line.split())>4:
            dip=[float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]

        #if read_orbital_energies and "-----" not in line and "..." not in line and len(line.split())!=0 and "Occupation" not in line:
            #num_entries=len(line.split())
            #if "(HOMO)" in line or "(LUMO)" in line:
            #    num_entries-=1
            #if not occ_done:
            #    if num_entries==4:
            #        occ_energies.append(float(line.split()[3]))
            #    elif num_entries==3:
            #        print("WARNING: error in parsing orbital energies")
            #        occ_energies=[]
            #        virt_energies=[]
            #        read_orbital_energies=False
            #else:
            #    if num_entries==4:
            #        if "0.0000" in line:
            #            pass
            #        else:
            #            print("WARNING: unexpected number of columns in parsing virtual energies")
            #            print(line)
            #        virt_energies.append(float(line.split()[3]))
            #    elif num_entries==3:
            #        virt_energies.append(float(line.split()[2]))
            #if "(HOMO)" in line:
            #    occ_done=True

        if "molecular dipole:" in line:
            print ('dip_line={}'.format(line))
            read_dip=True
        if ":: HOMO-LUMO gap" in line:
            HOMO_LUMO_gap=float(line.split()[3])
        if "Orbital Energies and Occupations" in line:
            read_orbital_energies=True
            occ_energies=[]
            virt_energies=[]
            occ_done=False

    if len(occ_energies)==0:
        occ_energies=None
    else:
        occ_energies=occ_energies[1:]
    if len(virt_energies)==0:
        virt_energies=None
    else:
        virt_energies=virt_energies[:-1]

    alphas_dict = {"metal_idx": None}
    covCN_dict = {"metal_idx": None}
    muls_dict = {"metal_idx": None}
    c6_dict = {"metal_idx": None}
                      
    if len(alphas)!=0:
        alphas_dict = {"metal": alphas[metal_idx], "br": alphas[br_idx]}
        covCN_dict = {"metal": covCN[metal_idx], "br": covCN[br_idx]}
        muls_dict = {"metal": muls[metal_idx],  "br": muls[br_idx]}
        c6_dict = {"metal": c6[metal_idx], "br": c6[br_idx]}

    return(True, muls_dict, alphas_dict, covCN_dict, c6_dict, wbo_values, dip, HOMO_LUMO_gap) #occ_energies, virt_energies)


def xtb_opt_save(moldir, coords, elements, charge=0, spin=0, opt='yes', freeze=[]): #YT added

    print ('opt={}'.format(opt))
    startdir=os.getcwd()
    os.chdir(moldir)
    exportXYZ(coords, elements, "in.xyz")

    if len(freeze)>0:
        outfile=open("xcontrol","w")
        outfile.write("$fix\n")
        outfile.write(" atoms: ")
        for counter,i in enumerate(freeze):
            if (counter+1)<len(freeze):
                outfile.write("%i,"%(i+1))
            else:
                outfile.write("%i\n"%(i+1))
        #outfile.write("$gbsa\n solvent=toluene\n")
        outfile.close()
        add=" -I xcontrol "
    else:
        add=""

    if opt == 'yes':
        if charge==0 and spin==0:
            os.system("xtb %s in.xyz --opt >> xtb.log"%(add))
        else:
            os.system("xtb %s in.xyz --opt --chrg %i --uhf %i >> xtb.log"%(add,charge,spin))
        if not os.path.exists("xtbopt.xyz"):
            print("WARNING: xtb geometry optimization did not work %s"%(moldir))
            os.chdir(startdir)
            os.system("rm -r %s"%(rundir))
            return(coords, elements)

    else:
        print ('no opt')
        if charge==0 and spin==0:
            os.system("xtb %s in.xyz >> xtb.log"%(add))
        else:
            os.system("xtb %s in.xyz --chrg %i --uhf %i >> xtb.log"%(add,charge,spin))
        shutil.copy('in.xyz', 'xtbopt.xyz')

    coords_new, elements_new=readXYZ("xtbopt.xyz")
    os.chdir(startdir)
    return(coords_new, elements_new)


def readXYZ(filename):
    infile=open(filename,"r")
    coords=[]
    elements=[]
    lines=infile.readlines()
    if len(lines)<3:
        exit("ERROR: no coordinates found in %s/%s"%(os.getcwd(), filename))
    for line in lines[2:]:
        elements.append(line.split()[0].capitalize())
        coords.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    infile.close()
    coords=np.array(coords)
    return coords,elements

def exportXYZ(coords,elements,filename, mask=[]):
    outfile=open(filename,"w")

    if len(mask)==0:
        outfile.write("%i\n\n"%(len(elements)))
        for atomidx,atom in enumerate(coords):
            outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(),atom[0],atom[1],atom[2]))
    else:
        outfile.write("%i\n\n"%(len(mask)))
        for atomidx in mask:
            atom = coords[atomidx]
            outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(),atom[0],atom[1],atom[2]))
    outfile.close()

def try_mkdir(dirname):
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            pass

def run_morfeus(coords, elements, moldir, settings):
    
    outfilename="%s/morfeus.yml"%(moldir)
    if os.path.exists(outfilename):
        infile=open(outfilename,"r")
        results = yaml.load(infile, Loader=yaml.FullLoader)
        infile.close()
        return(results)

    if settings['run_morfeus'] is False:
        return ('No metal, do not run morfeus')
    
    else:
        metal_idx = settings['metal_idx']+1 #morfeus is 1 based (metal_idx is 0 based)
        do_pyramid=True
        metal_char=settings["metal"]

        time_start=time.time()

        #Calculate cone_angle
        try:
            cone_angle = ConeAngle(elements, coords, metal_idx)
            cone_angle_val = float(cone_angle.cone_angle)
            print ('cone_angle={}'.format(cone_angle))
            print ('cone_angle_val={}'.format(cone_angle_val))
        except:
            print("WARNING: morfeus cone angle failed")
            cone_angle_val = None

        #Calculate solid_angle
        try:
            solid_angle = SolidAngle(elements, coords, metal_idx)
            solid_angle_val = float(solid_angle.solid_angle)
            solid_cone_angle_val = float(solid_angle.cone_angle)
            G_val = float(solid_angle.G)
            print ('solid_angle={}'.format(solid_angle))
            print ('solid_angle_val={}'.format(solid_angle_val))
        except:
            print("WARNING: morfeus solid angle failed")
            solid_angle_val = None
            solid_cone_angle_val = None
            G_val = None
            
        #Calculate SASA
        try:
            sasa = SASA(elements, coords)
            print ('sasa = {}'.format(sasa))
            sasa_val = float(sasa.area)
            sasa_val_metal = float(sasa.atom_areas[metal_idx])
            sasa_volume = float(sasa.volume)
            sasa_volume_metal = float(sasa.atom_volumes[metal_idx])
        except:
            print("WARNING: morfeus sasa failed")
            sasa_val = None
            sasa_val_metal = None
            sasa_volume = None
            sasa_volume_metal = None

        # Calculate Sterimol
        #try:
            #sterimol = Sterimol(elements_extended, coords_extended, dummy_idx+1, p_idx+1)
            #lval = float(sterimol.L_value)
            #B1 = float(sterimol.B_1_value)
            #B5 = float(sterimol.B_5_value)
        #except:
            #print("WARNING: morfeus sterimol failed")
            #lval = None
            #B1 = None
            #B5 = None

        # Calculate Dispersion, check if surface area is needed?
#        try:
        disp = Dispersion(elements, coords)
        print (disp.print_report())
        p_int = float(disp.p_int) #p_int: dispersion descriptor, based on electron density
        p_int_atoms_metal = float(disp.atom_p_int[metal_idx])
        p_int_area = float(disp.area)
        p_int_atom_areas_metal = float(disp.atom_areas[metal_idx])
        p_int_volume = float(disp.volume)
        print ('p_int={}'.format(p_int))
        print ('p_int_atoms_metal={}'.format(p_int_atoms_metal))
        print ('p_int_area={}'.format(p_int_area))
        print ('p_int_atom_areas_metal={}'.format(p_int_atom_areas_metal))
        print ('p_int_volume={}'.format(p_int_volume))

        #except:
        #    print("WARNING: morfeus dispersion failed")
        #    p_int = None
        #    p_int_atoms_metal = None
        #    p_int_area = None
        #    p_int_atom_areas_metal = None

        # Calculate Pyramidalization - two equivalent measurments P and alpha
        #if do_pyramid:
        #    try:
        #        pyr = Pyramidalization(elements = elements_extended, coordinates = coords_extended, atom_index = p_idx+1, excluded_atoms = [dummy_idx+1]) # remove Pd
        #        pyr_val = float(pyr.P)
        #        pyr_alpha = float(pyr.alpha)
        #    except:
        #        print("WARNING: morfeus Pyramidalization failed")
        #        pyr_val = None
         #       pyr_alpha = None
        #else:
        #    pyr_val = None
        #    pyr_alpha = None

        # Calculate Buried volume
        try:
            bv = BuriedVolume(elements, coords, metal_idx, density=0.01) # dummy_idx+1 = 2
            vbur = float(bv.buried_volume)
            #vdist = float(bv.distal_volume)  # the volume of the ligand beyond the sphere
            #vtot = float(vbur + vdist)
            vbur_percent = float(bv.fraction_buried_volume)
            print ('vbur={}'.format(vbur))
            print ('vbur_percent={}'.format(vbur_percent))
        except:
            print("WARNING: morfeus BuriedVolume failed")
            vbur = None
            #vtot = None
            vbur_percent = None

    time_end = time.time()
    time_morfeus = time_end - time_start
    
    results={#"lval": lval,
             #"B1": B1,
             #"B5": B5,
             "sasa": sasa_val,
             "sasa_metal": sasa_val_metal,
             "sasa_volume": sasa_volume,
             "sasa_volume_metal": sasa_volume_metal,
             "cone_angle": cone_angle_val,
             "solid_angle": solid_angle_val,
             "solid_cone_angle": solid_cone_angle_val,
             "solid_G_val": G_val,
             "p_int": p_int,
             "p_int_atoms_metal": p_int_atoms_metal,
             "p_int_area": p_int_area,
             "p_int_atom_areas_metal": p_int_atom_areas_metal,
             "p_int_volume": p_int_volume,
             #"pyr_val": pyr_val,
             #"pyr_alpha": pyr_alpha,
             "vbur": vbur,
             #"vtot": vtot,
             "vbur_percent": vbur_percent,
             "metal_idx": int(metal_idx),
             "time_morfeus": time_morfeus
             }

    outfilename="%s/morfeus.yml"%(moldir) 
    outfile=open(outfilename, "w")
    outfile.write(yaml.dump(results, default_flow_style=False))
    outfile.close()

    print (results)
    return(results)

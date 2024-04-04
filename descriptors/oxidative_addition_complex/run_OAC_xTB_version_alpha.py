# modified from krakenA
import os
import sys
import numpy as np
import subprocess
import shlex
import time
import yaml
import copy
from rdkit import Chem
# imported library
import oac_utils

def main(argv):

    startdir = os.getcwd()
    subdir = sys.argv[14]
    inputdir = "input_structures"+'/'+subdir
    if os.path.isdir(inputdir) is False:
        os.makedirs(inputdir)
    resultsdir = "results_all"+'/'+subdir
    if os.path.isdir(resultsdir) is False:
        os.makedirs(resultsdir)
    debugdir = "debugging"+'/'+subdir
    if os.path.exists(debugdir) is False:
        os.makedirs(debugdir)
    calcdir = "calculations"+'/'+subdir
    if os.path.exists(calcdir) is False:
        os.makedirs(calcdir)
        
    molname = sys.argv[2]
    chrg = int(sys.argv[6])
    spin = int(sys.argv[8])
    ligand_type = str(sys.argv[10])
    opt = sys.argv[12]

    print ('ligand_type={}'.format(ligand_type))

    if ligand_type == 'ferrocene-phosphine': #for ferrocene-phosphine, need to use constraint xtb, case by case 
        opt = 'no'
        
    print ('molname,chrg,spin={} {} {}'.format(molname, chrg, spin))
    coords, elements = oac_utils.readXYZ("%s/%s.xyz"%(inputdir, molname))
    os.chdir(calcdir)
    moldir = molname
    oac_utils.try_mkdir(moldir)

    time_start = time.time()
    # Check if already xtb done
    if os.path.exists(moldir+'/xtbopt.xyz'):
        print ('xtb already done')
        coords, elements = oac_utils.readXYZ(moldir+"/xtbopt.xyz")
    else:
        print ('run xtb')
        coords, elements = oac_utils.xtb_opt_save(moldir, coords, elements, charge=chrg, spin=spin, opt=opt, freeze=[])

    #coords_masked, elements_masked = 
        
    settings = {"metal": False,
                "run_morfeus": False
               }
    metals = ['Pd','Ni','Cu']
    for metal in metals:
        if metal in elements:
            settings["metal"] = metal
            settings["run_morfeus"] = True
            settings["metal_idx"] = elements.index(metal)
            break
    print (settings)

    settings["ligand"]=ligand_type
    adjmat = oac_utils.Table_generator(elements,coords,File=None)
    settings["metal_bond_count"] = int(sum(adjmat[settings["metal_idx"]]==1)) #check if bite_angle is necessary
    metal_bonded_atomids = np.where(adjmat[settings["metal_idx"]] == 1)[0]
    metal_bonded_elements = [elements[idx] for idx in metal_bonded_atomids]
    metal_bonded_map = {}
    for elem, atomid in zip(metal_bonded_elements, metal_bonded_atomids):
        if elem in metal_bonded_map:
            metal_bonded_map[elem].append(atomid)
        else:
            metal_bonded_map[elem] = [atomid]
    #print (metal_bonded_map)

    metal_idx = elements.index(metal)
    #br_idx = metal_bonded_map.get('Br', None)[0] #if non-bonding??
    br_idx = metal_bonded_map.get('Br', [None])[0]
    if br_idx is None:
        br_idx = next((i for i, element in enumerate(elements) if element == 'Br'), None)

    c_idx_list = metal_bonded_map.get('C', None)
    cf3ph_C_idx = None
    ligand_indices = []
    
    #Check NHC first
    if settings['ligand'] =='NHC':
        for idx in c_idx_list:
            c_bonded_atoms = np.where(adjmat[idx] == 1)[0]
            if sum(adjmat[idx]) == 3:
                c_bonded_elements = [elements[idx] for idx in c_bonded_atoms]
                if c_bonded_elements.count('N') == 2 and metal_idx in c_bonded_atoms:
                    ligand_indices = [idx]
                else:
                    cf3ph_C_idx = idx

    else:
        if len(c_idx_list) == 1:
            cf3ph_C_idx = c_idx_list[0]
        else:
            for idx in c_idx_list:
                if sum(adjmat[idx]) == 3:
                    cf3ph_C_idx = idx
    
    #Check other ligand
    if settings["ligand"] == 'monophosphine' or settings["ligand"] == 'ferrocene-phosphine' or settings['ligand'] == 'diphosphine':
        ligand_indices.extend(metal_bonded_map['P'])
    elif settings["ligand"] == 'monodentate-nitrogen' or settings["ligand"] == 'bidentate-nitrogen' or settings["ligand"] == 'multidentate-nitrogen':
        ligand_indices.extend(metal_bonded_map['N'])
    elif settings["ligand"] == 'dioxo':
        ligand_indices.extend(metal_bonded_map['O'])
    elif settings["ligand"] == 'oxo':
        ligand_indices.extend(metal_bonded_map['O'])
    elif settings['ligand'] == 'N,O-bidentate':
        if 'N' in metal_bonded_map:
            ligand_indices.extend(metal_bonded_map['N'])
        if 'O' in metal_bonded_map:
            ligand_indices.extend(metal_bonded_map['O'])
    else:
        ligand_indices = ligand_indices

    print (cf3ph_C_idx)
    print (ligand_indices)
    
    # Break the bond between cf3ph_C_idx and metal_idx
    adjmat[cf3ph_C_idx][settings['metal_idx']] = 0
    adjmat[settings['metal_idx']][cf3ph_C_idx] = 0

    # Identify fragments
    fragment1 = oac_utils.DFS(adjmat, cf3ph_C_idx)
    fragment2 = oac_utils.DFS(adjmat, settings['metal_idx'])
    coords_frag2 = [coords[i] for i in fragment2 if i!=br_idx]
    elements_frag2 = [elements[i] for i in fragment2 if i!=br_idx]
    oac_utils.exportXYZ(coords_frag2, elements_frag2,moldir+"/frag_metal-L.xyz")
    adjmat_frag2 = oac_utils.Table_generator(elements_frag2,coords_frag2,File=None)

    settings_frag2 = {"metal_idx": elements_frag2.index(metal),
                      "run_morfeus": settings["run_morfeus"],
                      "metal": settings["metal"]
                     }
    # run morfeus of xTB
    morfeus_parameters = oac_utils.run_morfeus(coords_frag2, elements_frag2, moldir, settings_frag2) ##xTB, only one morfeus_parameters
    
    idx_dict = {'metal_idx': metal_idx,
                'ligand_indices': ligand_indices,
                'br_idx': br_idx
                }
    
    # parse electronic properties of xtb
    xtb_done, muls, alphas, covCN, c6, wbo, dip, HOMO_LUMO_gap = oac_utils.read_xtb_log1(moldir,len(elements),idx_dict)
    electronic_properties = []
    electronic_properties.append({"muls":muls,
                                  "alphas": alphas,
                                  "covCN":covCN,
                                  "c6": c6,
                                  "dip":dip,
                                  "HOMO_LUMO_gap":HOMO_LUMO_gap,
                                  #"occ_energies":occ_energies,
                                  #"virt_energies":virt_energies,
                                  "wbo": wbo
                                  })

    print (electronic_properties)
    
    # go back to the start directory
    os.chdir(startdir)

    # save data
    print("   ---   Save the results of molecule %s to %s/%s.yml"%(molname,resultsdir, molname))
    if xtb_done: #crest_done and xtb_done:
        # start putting data in the results dictionary
        results_here={}
        results_here["coords"]=coords.tolist()
        results_here["elements"]=elements

        results_here["xTB"] = {"coords":[[float(val) for val in sublist] for sublist in coords], #float in case yaml dump error
                               "elements":elements,
                               "electronic_properties":electronic_properties,
                               "morfeus_parameters":morfeus_parameters  ###
                               }

        # add the timings
        time_end=time.time()
        time_all=time_end-time_start
        results_here["settings"]=settings
        results_here["time_all"]=time_all

        # save the main output file (this will hopefully be the smallest file with the most important data
        outfilename="%s/%s/%s.yml"%(startdir,resultsdir, molname)
        outfile=open(outfilename,"w")
        outfile.write(yaml.dump(results_here, default_flow_style=False))
        outfile.close()

    else:
        print("   ---   molecule %s FAILED"%(molname))
        outfilename="%s/%s/%s.yml"%(startdir,resultsdir, molname)
        outfile=open(outfilename,"w")
        outfile.write("FAILED\n")
        outfile.close()
        

    print("   ###   Finished molecule %s"%(molname))


if __name__ == "__main__":
    main(sys.argv[1:])

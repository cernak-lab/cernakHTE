# modified from kraken
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
import ligand_utils

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
    coords, elements = ligand_utils.readXYZ("%s/%s.xyz"%(inputdir, molname))
    os.chdir(calcdir)
    moldir = molname
    ligand_utils.try_mkdir(moldir)

    time_start = time.time()
    # Check if already xtb done
    if os.path.exists(moldir+'/xtbopt.xyz'):
        print ('xtb already done')
        coords, elements = ligand_utils.readXYZ(moldir+"/xtbopt.xyz")
    else:
        print ('run xtb')
        coords, elements = ligand_utils.xtb_opt_save(moldir, coords, elements, charge=chrg, spin=spin, opt=opt, freeze=[])

    settings = {"run_morfeus": False,
                "ligand": ligand_type
               }

    adjmat = ligand_utils.Table_generator(elements,coords,File=None)
    ligand_indices = []
    
    #Check NHC first
    if settings['ligand'] =='NHC':
        C_indices = ligand_utils.get_indices_for_element('C', elements)
        ligand_indices = [idx for idx in C_indices if len(np.where(adjmat[idx] == 1)[0]) == 2 and all(elements[j] == 'N' for j in np.where(adjmat[idx] == 1)[0])]
    #Check other ligands
    if settings["ligand"] == 'monophosphine' or settings["ligand"] == 'ferrocene-phosphine' or settings['ligand'] == 'diphosphine':
        ligand_indices=ligand_utils.get_indices_for_element('P', elements)
    if settings["ligand"] == 'monodentate-nitrogen' or settings["ligand"] == 'bidentate-nitrogen' or settings["ligand"] == 'multidentate-nitrogen':
        N_indices=ligand_utils.get_indices_for_element('N', elements)
        CN_indices = [idx for idx in N_indices if len(np.where(adjmat[idx] == 1)[0]) == 1 and 'C' in [elements[j] for j in np.where(adjmat[idx] == 1)[0]]]
        ligand_indices = [idx for idx in N_indices if idx not in CN_indices]
    if settings["ligand"] == 'dioxo' or settings["ligand"] == 'oxo':
        O_indices = ligand_utils.get_indices_for_element('O', elements)
        ligand_indices=[idx for idx in O_indices if len(np.where(adjmat[idx] == 1)[0]) == 1 and 'C' in [elements[j] for j in np.where(adjmat[idx] == 1)[0]]] #ketone identifier
    if settings['ligand'] == 'N,O-bidentate':
        N_indices=ligand_utils.get_indices_for_element('N', elements)
        print (N_indices)
        NH_indices = [idx for idx in N_indices if 'H' in [elements[j] for j in np.where(adjmat[idx] == 1)[0]]]
        print (NH_indices)
        O_indices = ligand_utils.get_indices_for_element('O', elements)
        print (O_indices)
        CO_indices = [idx for idx in O_indices if len(np.where(adjmat[idx] == 1)[0]) == 1 and 'C' in [elements[j] for j in np.where(adjmat[idx] == 1)[0]]]
        print (CO_indices)
        ligand_indices = NH_indices + CO_indices
        print (ligand_indices)
        
        
    print (ligand_indices)

    # parse ligand electronic properties of xtb
    #xtb_done, muls, alphas, covCN, c6, wbo, dip, HOMO_LUMO_gap = ligand_utils.read_xtb_log1(moldir,len(elements),idx_dict)
    xtb_done, muls, alphas, covCN, c6, wbo, dip, HOMO_LUMO_gap, E_HOMO, E_LUMO = ligand_utils.read_xtb_log2(moldir,len(elements),ligand_indices)
    electronic_properties = []
    electronic_properties.append({"muls":muls,
                                  "alphas": alphas,
                                  "covCN":covCN,
                                  "c6": c6,
                                  "dip":dip,
                                  "HOMO_LUMO_gap":HOMO_LUMO_gap,
                                  "E_HOMO":E_HOMO,
                                  "E_LUMO":E_LUMO,
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
                               #"morfeus_parameters":morfeus_parameters  ###
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

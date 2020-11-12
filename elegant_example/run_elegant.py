import sys
#lume_path = ''
#sys.path.append(lume_path + 'lume-impact')
#sys.path.append(lume_path + 'openPMD-beamphysics')
#sys.path.append(lume_path + 'distgen')
#sys.path.append(lume_path + 'xopt')


import numpy as np
import subprocess
from pmd_beamphysics.interfaces.elegant import elegant_h5_to_data
from pmd_beamphysics.plot import marginal_plot
from pmd_beamphysics import ParticleGroup
from pmd_beamphysics.plot import slice_statistics
import tempfile
import json
import shutil
import os
import h5py
from h5py import File
import re
from xopt.tools import fingerprint


def execute_elegant(settings, ele_fname = None, lte_fname = None):
    
    ####---------#### grab other relevant info for run from settings
    ####---------#### pop ones that aren't used at the moment too so are not in settings1 -- can fix later
    settings1 = settings.copy()
    ele_fname = settings1.pop('ele_fname')
    lte_fname = settings1.pop('lte_fname')
    path_search = settings1.pop('path_search')
    ELEGANT_BIN = settings1.pop('ELEGANT_BIN')
    HDF5_BIN = settings1.pop('HDF5_BIN')

    finput_name = settings1.pop('finput_name')
    foutput_name = settings1.pop('foutput_name')
    

    ####---------#### temp dir for run
    tdir = tempfile.TemporaryDirectory()
    path = tdir.name

    shutil.copy(ele_fname, path)
    shutil.copy(lte_fname, path)
    
    ####---------#### adjust search path according to YAML and rewrite .ele as 'run<ele_fname>' 
    
    with open(path+'/'+ele_fname,'r') as f:
    
        outlines=[]

        for line in f.readlines():

            x = line.strip()

            if x.startswith('search_path'):
                line = ' search_path = "'+path_search+'"\n'

            outlines.append(line)
            
    
    with open(path+'/run'+ele_fname, 'w') as f:
        for line in outlines:
            f.write(line)
    
    ####---------#### run elegant

    binname = ELEGANT_BIN
    cmd = binname + ' ' + 'run' + ele_fname 
    
    for s in settings1:
        cmd += ' -macro=' + s + '=' + str(settings1[s]) + ' '
        
    cmd_run = cmd.split()
    output = {'error':True, 'log':''}

    p = subprocess.run(cmd_run, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=path)
    
    ####---------#### process output
    PF = process_output(path,HDF5_BIN,finput_name,foutput_name)
    output['log'] = p.stdout
    output['error'] = False

    return PF

def process_output(path, HDF5_BIN, finput_name = None, foutput_name = None,timeout=None):

    finput = os.path.join(path,finput_name)
    foutput = os.path.join(path,foutput_name)

    cmd = f'{HDF5_BIN} {finput} {foutput}'
    cmd_run = cmd.split()
    output = {'error':True, 'log':''}

    p = subprocess.run(cmd_run, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    data = elegant_h5_to_data(foutput)
    PF = ParticleGroup(data=data)

    output['log'] = p.stdout
    output['error'] = False
    
    return PF


def calculate_core_pars(P, n_slice=50, n_slice_av=10):
    
    n_lower = int((n_slice - n_slice_av)/2)
    n_upper = int((n_slice + n_slice_av)/2)
        
    P.t -= P['mean_t']
    keys = ['norm_emit_x', 'sigma_pz', 'ptp_t','charge']
    slice_dat = slice_statistics(P, keys, n_slice, 't')
    slice_dat['density'] = slice_dat['charge']/ slice_dat['ptp_t']
    
    core_current = np.mean(slice_dat['density'][n_lower:n_upper])
    core_pz = np.mean(slice_dat['sigma_pz'][n_lower:n_upper])
    core_norm_emit_x = np.mean(slice_dat['norm_emit_x'][n_lower:n_upper])

    core_dict = {'core_norm_emit_x':core_norm_emit_x, 'core_pz':core_pz, 'core_current': core_current}

    return core_dict


def merit1(P):

    d = {}
    for k in ['sigma_pz', 'norm_emit_x', 'sigma_t', 'mean_pz', 'charge', 'ptp_t']:
        d[k] = P[k]
     
    core_dict = calculate_core_pars(P)
    d_up = {**d, **core_dict}
    
    return d_up

def evaluate_elegant(settings, a=1):
    
    H5_SAVE = settings.pop('H5_SAVE')
    
    ####---------#### run elegant
    
    PF = execute_elegant(settings)
    
    ####---------#### save h5 and archive
    
    fnameh5 = 'elegant_sim_' + fingerprint(settings) + '.h5'
    PF.write(H5_SAVE+fnameh5)

    output = merit1(PF)
    output['archive']=H5_SAVE+fnameh5
 
    return output


if __name__ == 'main':

    settings = {'L1_9_50_phase':50.5, 'L1_9_25_phase':50.5, 'L1_10_25_phase':50.5, 'X1_Xband_phase':-70, 
                'L2_10_50_phase':55.5, 'L2_10_25_phase':55.5, 'L3_10_25_volt':1.6628471874e7, 'X_MAX':1.3e-3,             
                'DX':0.9e-3, 'DP': 15.0e-5,
                'INPUT_FILE':'elegant_particles.txt'}
    
    PF = execute_elegant(settings)


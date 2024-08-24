
"""
netParams.py

High-level specifications for S1-thalamus network model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle, json
import os
import numpy as np
import pandas as pd

netParams = specs.NetParams()   # object of class NetParams to store the network parameters


try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------
netParams.scale = cfg.scale # Scale factor for number of cells
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume

# L23 Human net
# #              L2/3   L4     L5
PYRmaxApics = [550   ,1550   ,1900]
uppers =      [-250  ,-1200 ,-1600]
lowers =      [-1200 ,-1580 ,-2300]

L25_human = 250 + 950 + 380 + 720 + 1000
Human_height = 3300.0

cellModels = ['HH_full']

Ipops = ['HL23VIP', 'HL23PV', 'HL23SST']
Epops = ['HL23PYR']

layer = {'1':[0.0, 250.0], '23': [250.0,1200.0], '23soma': [550.0,1200.0], '4':[1200.0,1580.0], '5': [1580.0,2300.0], '6': [2300.0,3300.0]}  # normalized layer boundaries

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 0.1 # default conn delay (ms)
netParams.propVelocity = 300.0 #  300 μm/ms (Stuart et al., 1997)
netParams.rotateCellsRandomly = True

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------

for cellName in cfg.allpops:
    cellRule = netParams.importCellParams(label=cellName, somaAtOrigin=True,
        conds={'cellType': cellName, 'cellModel': 'HH_full'},
        fileName='cellwrapper.py',
        cellName='loadCell_' + cellName,
        cellInstance = True,
        cellArgs={'cellName': cellName})
    
#------------------------------------------------------------------------------
#Import Excel file
circuit_params = pd.read_excel('Circuit_param.xls', sheet_name = None, index_col = 0)

#Get cell names and import biophys
cell_names = [i for i in circuit_params['conn_probs'].axes[0]]

circuit_params["syn_params"] = {'none':{'tau_r_AMPA': 0,'tau_d_AMPA': 0,'tau_r_NMDA': 0,
                                'tau_d_NMDA': 0, 'e': 0,'Dep': 0,'Fac': 0,'Use': 0,'u0':0,'gmax': 0}}
circuit_params["multi_syns"] = {'none':{'loc':0,'scale':0}}
# organizing dictionary for LFPY input
for pre in cell_names:
    for post in cell_names:
        if "PYR" in pre:
            circuit_params["syn_params"][pre+post] = {'tau_r_AMPA': 0.3, 'tau_d_AMPA': 3, 'tau_r_NMDA': 2,
                                                      'tau_d_NMDA': 65, 'e': 0, 'u0':0,
                                                      'Dep': circuit_params["Depression"].at[pre, post],
                                                      'Fac': circuit_params["Facilitation"].at[pre, post],
                                                      'Use': circuit_params["Use"].at[pre, post],
                                                      'gmax': circuit_params["syn_cond"].at[pre, post]}
        else:
            circuit_params["syn_params"][pre+post] = {'tau_r': 1, 'tau_d': 10, 'e': -80, 'u0':0,
                                                      'Dep': circuit_params["Depression"].at[pre, post],
                                                      'Fac': circuit_params["Facilitation"].at[pre, post],
                                                      'Use': circuit_params["Use"].at[pre, post],
                                                      'gmax': circuit_params["syn_cond"].at[pre, post]}
        circuit_params["multi_syns"][pre+post] = {'loc':int(circuit_params["n_cont"].at[pre, post]),'scale':0}

stimuli = []
for stimulus in circuit_params['STIM_PARAM'].axes[0]:
    stimuli.append({})
    for param_name in circuit_params['STIM_PARAM'].axes[1]:
        stimuli[-1][param_name] = circuit_params['STIM_PARAM'].at[stimulus, param_name]
    new_param = circuit_params["syn_params"][stimuli[-1]['syn_params']].copy()
    new_param['gmax'] = stimuli[-1]['gmax']
    stimuli[-1]['new_param'] = new_param

#------------------------------------------------------------------------------
from scipy import stats as st

halfnorm_rv = st.halfnorm
uniform_rv = st.uniform

#              L2/3   L4     L5
PYRmaxApics = [550   ,1550   ,1900]
uppers =      [-250  ,-1200 ,-1600]
lowers =      [-1200 ,-1580 ,-2300]

depths = []
rangedepths = []
minSynLocs = []
syn_pos = []
pop_args = {}

for i in range (3):
    depths.append((lowers[i]-uppers[i])/2-PYRmaxApics[i])
    rangedepths.append(abs(lowers[i]-uppers[i])/2)
    minSynLocs.append((lowers[i]-uppers[i])/2*3-PYRmaxApics[i])

    syn_pos.append({'section' : ['apic', 'dend'],
                    'fun' : [uniform_rv, halfnorm_rv],
                    'funargs' : [{'loc':minSynLocs[i], 'scale':abs(minSynLocs[i])},{'loc':minSynLocs[i], 'scale':abs(minSynLocs[i])}],
                    'funweights' : [1, 1.]})
    syn_pos.append({'section' : ['apic'],
                    'fun' : [uniform_rv],
                    'funargs' : [{'loc':minSynLocs[i], 'scale':abs(minSynLocs[i])}],
                    'funweights' : [1.]})
    syn_pos.append({'section' : ['dend'],
                    'fun' : [uniform_rv],
                    'funargs' : [{'loc':minSynLocs[i], 'scale':abs(minSynLocs[i])}],
                    'funweights' : [1.]})
    syn_pos.append({'section' : ['dend'],
                   'fun' : [halfnorm_rv],
                   'funargs' : [{'loc':minSynLocs[i], 'scale':abs(minSynLocs[i])}],
                   'funweights' : [1.]})
    names = ['HL2','HL4','HL5']
    pop_args[names[i]]={'radius':250,
                        'loc':depths[i],
                        'scale':rangedepths[i]*4,
                        'cap':rangedepths[i]}
    
#------------------------------------------------------------------------------

# Rotate to z as vertical
rotate_x = {}
rotate_y = {}
rotate_z = {}
rotate_x['HL23PYR'], rotate_x['HL23SST'], rotate_x['HL23PV'], rotate_x['HL23VIP'] = 1.57, 1.77, 1.26, -1.57
rotate_y['HL23PYR'], rotate_y['HL23SST'], rotate_y['HL23PV'], rotate_y['HL23VIP'] = 2.62, 2.77, 2.57, 3.57
rotate_z['HL23PYR'], rotate_z['HL23SST'], rotate_z['HL23PV'], rotate_z['HL23VIP'] = 0.0, 0.0, 0.0, 0.0

for cellName in netParams.cellParams.keys():

    cellType = netParams.cellParams[cellName]['conds']['cellType']

    x = rotate_x[cellType]
    y = rotate_y[cellType]
    z = rotate_z[cellType]

    for sectName in netParams.cellParams[cellName]['secs'].keys():

        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        theta = -x
        rotation_x = np.array([[1, 0, 0],
                                       [0, np.cos(theta), -np.sin(theta)],
                                       [0, np.sin(theta), np.cos(theta)]])
        
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_x)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams


        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        phi = -y
        rotation_y = np.array([[np.cos(phi), 0, np.sin(phi)],
                                       [0, 1, 0],
                                       [-np.sin(phi), 0, np.cos(phi)]])
        
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_y)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams


        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        gamma = -z
        rotation_z = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                                       [np.sin(gamma), np.cos(gamma), 0],
                                       [0, 0, 1]])
    
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_z)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams

# Rotate to Y as vertical axis
rotate_x = {}
rotate_y = {}
rotate_z = {}
rotate_x['HL23PYR'], rotate_x['HL23SST'], rotate_x['HL23PV'], rotate_x['HL23VIP'] = -1.5708, -1.5708, -1.5708, -1.5708
rotate_y['HL23PYR'], rotate_y['HL23SST'], rotate_y['HL23PV'], rotate_y['HL23VIP'] = 0.0, 0.0, 0.0, 0.0
rotate_z['HL23PYR'], rotate_z['HL23SST'], rotate_z['HL23PV'], rotate_z['HL23VIP'] = 0.0, 0.0, 0.0, 0.0

for cellName in netParams.cellParams.keys():

    cellType = netParams.cellParams[cellName]['conds']['cellType']

    x = rotate_x[cellType]
    y = rotate_y[cellType]
    z = rotate_z[cellType]

    for sectName in netParams.cellParams[cellName]['secs'].keys():

        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        theta = -x
        rotation_x = np.array([[1, 0, 0],
                                       [0, np.cos(theta), -np.sin(theta)],
                                       [0, np.sin(theta), np.cos(theta)]])
        
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_x)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams

        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        phi = -y
        rotation_y = np.array([[np.cos(phi), 0, np.sin(phi)],
                                       [0, 1, 0],
                                       [-np.sin(phi), 0, np.cos(phi)]])
        
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_y)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams

        sectParams_new = netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d']
        sectParams = []

        gamma = -z
        rotation_z = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                                       [np.sin(gamma), np.cos(gamma), 0],
                                       [0, 0, 1]])
    
        for i in range(len(sectParams_new)):
            x3d, y3d, z3d, L3d = sectParams_new[i]
            rel_pos = x3d, y3d, z3d

            # print(rel_pos)        
            rel_pos = np.dot(rel_pos, rotation_z)
            # print(rel_pos)
            pt3d = (rel_pos[0],rel_pos[1] , rel_pos[2], L3d)
            sectParams.append(pt3d)

        netParams.cellParams[cellName]['secs'][sectName]['geom']['pt3d'] = sectParams

# Change axon names
for cellName in netParams.cellParams.keys():
        
        axon_pt3d_x, axon_pt3d_y, axon_pt3d_z, soma_pt3d_diam =  netParams.cellParams[cellName]['secs']['soma_0']['geom']['pt3d'][-1]
        axon_pt3d_diam =  netParams.cellParams[cellName]['secs']['axon_0']['geom']['diam']
        axon_pt3d_L =  netParams.cellParams[cellName]['secs']['axon_0']['geom']['L']

        netParams.cellParams[cellName]['secs']['axon_0']['geom']['pt3d'] = [(axon_pt3d_x, axon_pt3d_y, axon_pt3d_z, axon_pt3d_diam),
                                                                          (axon_pt3d_x, axon_pt3d_y+axon_pt3d_L/2.0, axon_pt3d_z, axon_pt3d_diam),
                                                                          (axon_pt3d_x, axon_pt3d_y+axon_pt3d_L, axon_pt3d_z, axon_pt3d_diam)]

        axon1_pt3d_x, axon1_pt3d_y, axon1_pt3d_z, soma_pt3d_diam =  netParams.cellParams[cellName]['secs']['axon_0']['geom']['pt3d'][-1]
        axon1_pt3d_diam =  netParams.cellParams[cellName]['secs']['axon_1']['geom']['diam']
        axon1_pt3d_L =  netParams.cellParams[cellName]['secs']['axon_1']['geom']['L']

        netParams.cellParams[cellName]['secs']['axon_1']['geom']['pt3d'] = [(axon1_pt3d_x, axon1_pt3d_y, axon1_pt3d_z, axon1_pt3d_diam),
                                                                          (axon1_pt3d_x, axon1_pt3d_y+axon1_pt3d_L/2.0, axon1_pt3d_z, axon1_pt3d_diam),
                                                                          (axon1_pt3d_x, axon1_pt3d_y+axon1_pt3d_L, axon1_pt3d_z, axon1_pt3d_diam)] 
        
        if 'myelin_0' in netParams.cellParams[cellName]['secs'].keys():

                if 'myelin_0' not in netParams.cellParams[cellName]['secLists']['all']:
                        netParams.cellParams[cellName]['secLists']['all'].append('myelin_0')

                if 'myelin_0' not in netParams.cellParams[cellName]['secLists']['axonal']:                        
                        netParams.cellParams[cellName]['secLists']['axonal'].append('myelin_0')

                myelin0_pt3d_x, myelin0_pt3d_y, myelin0_pt3d_z, soma_pt3d_diam =  netParams.cellParams[cellName]['secs']['axon_1']['geom']['pt3d'][-1]
                myelin0_pt3d_diam =  netParams.cellParams[cellName]['secs']['myelin_0']['geom']['diam']
                myelin0_pt3d_L =  netParams.cellParams[cellName]['secs']['myelin_0']['geom']['L']

                netParams.cellParams[cellName]['secs']['myelin_0']['geom']['pt3d'] = [(myelin0_pt3d_x, myelin0_pt3d_y, myelin0_pt3d_z, myelin0_pt3d_diam),
                                                                                (myelin0_pt3d_x, myelin0_pt3d_y+myelin0_pt3d_L/2.0, myelin0_pt3d_z, myelin0_pt3d_diam),
                                                                                (myelin0_pt3d_x, myelin0_pt3d_y+myelin0_pt3d_L, myelin0_pt3d_z, myelin0_pt3d_diam)] 
# print and rename
for cellName in netParams.cellParams.keys():
    if 'myelin_0' in netParams.cellParams[cellName]['secs'].keys():

        netParams.renameCellParamsSec(label=cellName, oldSec='myelin_0', newSec='axon_2')      
            
        for secname2 in netParams.cellParams[cellName]['secLists'].keys():
            if 'myelin_0' in netParams.cellParams[cellName]['secLists'][secname2]:
                print('old ->',cellName,secname2,netParams.cellParams[cellName]['secLists'][secname2][-1])
                netParams.cellParams[cellName]['secLists'][secname2][-1] = 'axon_2'    
                print('new ->',cellName,secname2,netParams.cellParams[cellName]['secLists'][secname2][-1])

# print and include 'secLists' -> 'spiny'
for cellName in netParams.cellParams.keys():
    print(netParams.cellParams[cellName]['secLists'].keys())
    print(netParams.cellParams[cellName]['secLists']['basal'])
    print(netParams.cellParams[cellName]['secLists']['apical'])
    netParams.cellParams[cellName]['secLists']['spiny'] = {}
    nonSpiny = netParams.cellParams[cellName]['secLists']['axonal']
    nonSpiny.append('soma_0')
    netParams.cellParams[cellName]['secLists']['spiny'] = [sec for sec in netParams.cellParams[cellName]['secLists']['all'] if sec not in nonSpiny]
    print(netParams.cellParams[cellName]['secLists']['spiny'])

#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------
for cellName in cfg.allpops:
    netParams.popParams[cellName] = {'cellType': cellName, 'cellModel': 'HH_full', 'numCells': cfg.cellNumber[cellName], 'yRange': layer['23soma']} 

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
     for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}

        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}

        
# #------------------------------------------------------------------------------
# #  extracellular mechs
# #------------------------------------------------------------------------------
if cfg.addExternalStimulation:
    for celltyp in netParams.cellParams.keys():
        for secname in netParams.cellParams[celltyp]['secs'].keys():
            netParams.cellParams[celltyp]['secs'][secname]['mechs']['extracellular'] = {}
else:
    for celltyp in netParams.cellParams.keys():
        for secname in netParams.cellParams[celltyp]['secs'].keys():
            if 'extracellular' in netParams.cellParams[celltyp]['secs'][secname]['mechs'].keys():
                # print(celltyp, secname)
                del netParams.cellParams[celltyp]['secs'][secname]['mechs']['extracellular']


# ['HL23PYR', 'HL23SST', 'HL23PV', 'HL23VIP']

# Synaptic mechanism parameters
netParams.synMechParams['AMPA'] = {'mod': 'Exp2Syn', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}  # AMPA
netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 0.15, 'tau2': 15, 'e': 0}  # NMDA
netParams.synMechParams['GABAA'] = {'mod': 'Exp2Syn', 'tau1': 0.07, 'tau2': 9.1, 'e': -80}  # GABAA
netParams.synMechParams['GABAB'] = {'mod': 'Exp2Syn', 'tau1': 0.07, 'tau2': 9.1, 'e': -80}  # GABAB

weightStim = {}
weightStim['HL23PYR'] = 0.0002
weightStim['HL23SST'] = 0.0002
weightStim['HL23PV'] = 0.0005
weightStim['HL23VIP'] = 0.0004

rateStim = {}
rateStim['HL23PYR'] = 100.0
rateStim['HL23SST'] = 100.0
rateStim['HL23PV'] = 100.0
rateStim['HL23VIP'] = 100.0

SourcesNumber = 5 # for each post Mtype - sec distribution

for post in ['HL23PYR', 'HL23VIP', 'HL23PV', 'HL23SST']:

    for qSnum in range(SourcesNumber):
            
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (rateStim[post])/SourcesNumber

            netParams.stimSourceParams['bkg_' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': 100, 'noise': 1.0, 'start': 0}

            netParams.stimTargetParams['bkg->' + post + '_' + str(qSnum)] = {'source': 'bkg_' + post + '_' + str(qSnum), 'conds': {'pop': [post]}, 'weight': weightStim[post], 'delay': 0.5}


#------------------------------------------------------------------------------
# Connectivity parameters
#------------------------------------------------------------------------------

print(depths, minSynLocs,rangedepths,syn_pos[0]['section'],pop_args['HL2'],'\n\n')

for pre in cell_names:
    for post in cell_names:
        
        Syn_pos = int(circuit_params['Syn_pos'].at[pre, post])


        if circuit_params['conn_probs'].at[pre, post] > 0.0:

            if "PYR" in pre: # Excitatory
                netParams.synMechParams[pre+post] = {'mod': 'ProbAMPANMDA_EMS',
                                                        'tau_r_AMPA': 0.3, 'tau_d_AMPA': 3.0, 'tau_r_NMDA': 2.0,
                                                        'tau_d_NMDA': 65.0,
                                                        'Dep': circuit_params["Depression"].at[pre, post],
                                                        'Fac': circuit_params["Facilitation"].at[pre, post],
                                                        'Use': circuit_params["Use"].at[pre, post]}
            else:
                netParams.synMechParams[pre+post] = {'mod': 'ProbGABAAB_EMS',
                                                    'tau_r_GABAA': 1, 'tau_d_GABAA': 10,
                                                        'Dep': circuit_params["Depression"].at[pre, post],
                                                        'Fac': circuit_params["Facilitation"].at[pre, post],
                                                        'Use': circuit_params["Use"].at[pre, post]}
                

            netParams.connParams[pre + '->' + post] = {'preConds': {'cellType': pre}, 
                                                        'postConds': {'cellType': post},  #  E -> all (100-1000 um) ,'y': [0,5000]
                                                        'probability': circuit_params['conn_probs'].at[pre, post],                  # probability of connection
                                                        'weight': circuit_params['syn_cond'].at[pre, post],         # synaptic weight 
                                                        'delay': 0.5,      # transmission delay (ms) 
                                                        'synMech': pre+post,
                                                        'synsPerConn': int(circuit_params['n_cont'].at[pre, post]),
                                                        'sec': 'spiny',
                                                        }  
             
            if "PYR" in pre and "PYR" in post:
                netParams.connParams[pre + '->' + post]['sec'] = 'spiny'
            else:
                if 'dend' in syn_pos[Syn_pos]['section']:
                    netParams.connParams[pre + '->' + post]['sec'] = 'basal'
                if 'apic' in syn_pos[Syn_pos]['section']:
                    netParams.connParams[pre + '->' + post]['sec'] = 'apical'

            print('\n',pre+post,circuit_params['n_cont'].at[pre, post],circuit_params['conn_probs'].at[pre, post],Syn_pos,syn_pos[Syn_pos]['section'])
            print(netParams.synMechParams[pre+post])
            print(netParams.connParams[pre + '->' + post])
            print(netParams.connParams[pre + '->' + post]['sec'])

#------------------------------------------------------------------------------
# NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
#------------------------------------------------------------------------------
# SourcesNumber = 5 # for each post Mtype - sec distribution
# synperNeuronStimI = connData['synperNeuronStimI']
# synperNeuronStimE = connData['synperNeuronStimE']
# GsynStimI = connData['GsynStimI']
# GsynStimE = connData['GsynStimE']
   
# if cfg.addStimSynS1:      
#     for post in Ipops + Epops:

#         synperNeuron = synperNeuronStimI[post]
#         ratespontaneous = cfg.rateStimI
#         for qSnum in range(SourcesNumber):
#             ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
#             netParams.stimSourceParams['StimSynS1_S_all_INH->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}

#         synperNeuron = synperNeuronStimE[post]
#         ratespontaneous = cfg.rateStimE
#         for qSnum in range(SourcesNumber):
#             ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
#             netParams.stimSourceParams['StimSynS1_S_all_EXC->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}
            
#     #------------------------------------------------------------------------------
#     for post in Epops:
#         for qSnum in range(SourcesNumber):
#             netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
#                 'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
#                 'conds': {'cellType': cfg.popLabelEl[post]}, 
#                 'synMech': 'AMPA', 
#                 'sec': 'spinyEE', 
#                 'weight': GsynStimE[post],
#                 'delay': 0.1}

#     for post in Ipops:
#         for qSnum in range(SourcesNumber):
#             netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
#                 'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
#                 'synMech': 'AMPA', 
#                 'conds': {'cellType': cfg.popLabelEl[post]}, 
#                 'sec': 'spiny', 
#                 'weight': GsynStimE[post],
#                 'delay': 0.1}

#     for post in Epops+Ipops:
#         for qSnum in range(SourcesNumber):
#             netParams.stimTargetParams['StimSynS1_T_all_INH->' + post + '_' + str(qSnum)] = {
#                 'source': 'StimSynS1_S_all_INH->' + post + '_' + str(qSnum), 
#                 'conds': {'cellType': cfg.popLabelEl[post]}, 
#                 'synMech': 'GABAA', 
#                 'sec': 'spiny', 
#                 'weight': GsynStimI[post],
#                 'delay': 0.1}


#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
netParams.description = """ 
- Code based: LFPy version zenodo: ..... 
- v0 - L23 firing patterns
- v1 - include TMS and TACS

"""

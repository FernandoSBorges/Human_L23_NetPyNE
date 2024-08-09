"""
cfg.py 

Simulation configuration for L23 human model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: ... , fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import os
import numpy as np

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.simType='L23'
cfg.coreneuron = False

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 1000.0 ## Duration of the sim, in ms  
cfg.dt = 0.025
cfg.seeds = {'cell': 4321, 'conn': 4321, 'stim': 1234, 'loc': 4321} 
cfg.hParams = {'celsius': 34, 'v_init': -80}  
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True  
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1

cfg.includeParamsLabel = False
cfg.printPopAvgRates = True
cfg.checkErrors = False

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
cfg.reducedtest = False    

# TO DEBUG - Create only 5 Cells for each MEtype in S1
cfg.oneCellperMEtypeS1 = False 


cfg.rootFolder = os.getcwd()

# Load cells info from previously saved using netpyne (False: load from HOC BBP files, slower)
cfg.loadcellsfromJSON = True


cfg.poptypeNumber = 4
cfg.celltypeNumber = 4

cfg.allpops = ['HL23PYR', 'HL23SST', 'HL23PV', 'HL23VIP']

#--------------------------------------------------------------------------
# Recording 
#--------------------------------------------------------------------------
cfg.cellsrec = 0
if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop

cfg.recordTraces = {'V_soma': {'sec':'soma_0', 'loc':0.5, 'var':'v'},
                    'V_axon_0': {'sec':'axon_0', 'loc':0.5, 'var':'v'},                  
                    'V_dend_5': {'sec':'dend_5', 'loc':0.5, 'var':'v'},
                    }

cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.1           

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'v1_batch0'
cfg.saveFolder = 'Circuit_output/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = True         	## Save pkl file
cfg.saveJson = False	           	## Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net'] ## , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = True			
cfg.saveCellConns = True	

#------------------------------------------------------------------------------
# Analysis and plotting 
# ------------------------------------------------------------------------------

cfg.analysis['plotTraces'] = {'include': cfg.allpops, 'figSize': (12, 4), 'timeRange': [0,cfg.duration], 'saveFig': True, 'overlay': True, 'oneFigPer': 'cell', 'figSize':(12,4)}  # Plot recorded traces for this list of cells
cfg.analysis['plotRaster'] = {'include': cfg.allpops, 'saveFig': True, 'showFig': False, 'orderInverse': True, 'timeRange': [0,cfg.duration], 'figSize': (12,6), 
                              'fontSize':12, 'lw': 4, 'markerSize':4, 'marker': '.', 'dpi': 100} 
cfg.analysis['plot2Dnet']   = {'include': cfg.allpops, 'saveFig': True, 'showConns': True, 'figSize': (12,12), 'fontSize':12}   # Plot 2D cells xy
cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'oneFigPer': 'cell', 'overlay': True, 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False, 'figSize':(18,12)}
cfg.analysis['plotShape'] = {'includePre': cfg.recordCells, 'includePost': cfg.recordCells, 'showFig': False, 'includeAxon': False, 
                            'showSyns': False, 'saveFig': True, 'dist': 0.55, 'cvar': 'voltage', 'figSize': (24,24), 'dpi': 300}

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------

PYRmaxApics = [550   ,1550   ,1900]
uppers =      [-250  ,-1200 ,-1600]
lowers =      [-1200 ,-1580 ,-2300]

L25_human = 250 + 950 + 380 + 720 + 1000
Human_height = 3300.0

cfg.scale = 1.0 # reduce size
cfg.sizeY = 3300.0
cfg.sizeX = 250.0 # r =  um 
cfg.sizeZ = 250.0

cell_num = [800, 50, 70, 80]
# cell_num = [80, 5, 7, 8]

cfg.cellNumber = {}
for ii,cellName in enumerate(['HL23PYR', 'HL23SST', 'HL23PV', 'HL23VIP']):
    cfg.cellNumber[cellName] = cell_num[ii]

#------------------------------------------------------------------------------
# Spontaneous synapses + background
#------------------------------------------------------------------------------
cfg.addStimSynS1 = True
cfg.rateStimE = 100.0
cfg.rateStimI = 100.0

#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
## S1->S1
cfg.addConn = True

cfg.EEGain = 1.0
cfg.EIGain = 1.0
cfg.IIGain = 1.0
cfg.IEGain = 1.0

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 0

cfg.IClamp1 = {'pop': 'HL23PYR',  'sec': 'soma_0', 'loc': 0.5, 'start': 0, 'dur': 100, 'amp': 0.1}
# cfg.IClamp2 = {'pop': 'HL23VIP', 'sec': 'soma_0', 'loc': 0.5, 'start': 100, 'dur': 300, 'amp': 0.1}
# cfg.IClamp3 = {'pop': 'HL23PV', 'sec': 'soma_0', 'loc': 0.5, 'start': 100, 'dur': 300, 'amp': 0.1}
# cfg.IClamp4 = {'pop': 'HL23SST', 'sec': 'soma_0', 'loc': 0.5, 'start': 100, 'dur': 300, 'amp': 0.1}

#------------------------------------------------------------------------------
# External Stimulation
#------------------------------------------------------------------------------

cfg.addExternalStimulation = False

# The parameters of Alternate Current Stimulation
cfg.acs_params = {'position': [0.0, -1710.0, 0.0],  # um # y = [pia, bone]
              'amp': -1250.,  # uA,
              'stimstart': 300,  # ms
              'stimend': 400.0,  # ms
              'frequency': 5,  # Hz
              'sigma': 0.57  # decay constant S/m
              }

# The parameters of Transcranial Magnetic Stimulation 
cfg.tms_params = dict(
    freq_Hz=30.,
    duration_ms=cfg.duration,
    pulse_resolution_ms=cfg.dt,
    stim_start_ms=2000.,
    stim_end_ms=3000.,
    ef_amp_V_per_m=80.,
    width_ms=1.,
    pshape="Sine",
    decay_rate_percent_per_mm=10,
    E_field_dir=[-1, -1, -1],
    decay_dir=[0, 0, -1],
    ref_point_um=[0, 0, 0],
)

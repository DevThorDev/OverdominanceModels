# ### INPUT ####################################################################
import os, sys

ND_PROGFILES = 'ProgramFiles'
ND_TOOLS = 'Tools'
PD_PROGFILES = '../' + ND_PROGFILES
PD_TOOLS = '../' + ND_TOOLS
lDSub = [PD_PROGFILES, PD_TOOLS]
for cDSub in lDSub:
    sys.path.append(sys.path[0] + '/' + cDSub)
print('Current dir.:', os.getcwd())
import Constants as Cst

# flags
doCalcs = True
extractInfo = True
plotModelCmp = False

# additional model to compare against model defined by Cst.D_NRASOD_1
addModel = Cst.D_RASOD_1
# task list
# lDTasksIn = [{'minAFit': [0.0, 0.9],
#               'maxAFit': [0.1, 1.0],
#               'numAllIni': [50, 100, 250, 500, 1000],
#               'multNumPatho': [2, 10]},         # minimum: 2; even number
#              {'minAFit': [0.2],
#               'maxAFit': [0.6],
#               'numAllIni': [50, 100, 250, 500, 1000],
#               'multNumPatho': [2, 10]},         # minimum: 2; even number
#              {'minAFit': [0.4],
#               'maxAFit': [0.8],
#               'numAllIni': [50, 100, 250, 500, 1000],
#               'multNumPatho': [2, 10]},         # minimum: 2; even number
#              {'minAFit': [0.45],
#               'maxAFit': [0.55],
#               'numAllIni': [50, 100, 250, 500, 1000],
#               'multNumPatho': [2, 10]},         # minimum: 2; even number
#              {'minAFit': [0.3],
#               'maxAFit': [0.7],
#               'numAllIni': [50, 100, 250, 500, 1000],
#               'multNumPatho': [2, 10]}]         # minimum: 2; even number
# lDTasksIn = [{'minAFit': [0.0],
#               'maxAFit': [1.0],
#               'numAllIni': [5, 20],
#               'multNumPatho': [2]},         # minimum: 2; even number
#              {'minAFit': [0.3],
#               'maxAFit': [0.7],
#               'numAllIni': [50],
#               'multNumPatho': [2, 10]}]         # minimum: 2; even number
lDTasksIn = [{'minAFit': [0.0],
              'maxAFit': [1.0],
              'numAllIni': [100],
              'multNumPatho': [2]}]         # minimum: 2; even number
taskMode = Cst.N_TMD_MULT                   # N_TMD_MULT

# other input variables or constants
CALC_MODE = Cst.C_SOLVE_LE          # C_SOLVE_LE / C_TSTEPPING
IN_ALL_FDIST = Cst.A_FDIST_UNIFORM  # A_FDIST_UNIFORM / A_FDIST_RANDOM
IN_ALL_PROP = Cst.A_PROP_UNIFORM    # A_PROP_UNIFORM / A_PROP_RANDOM
DOM_TYPE = Cst.D_NRASOD_1           # D_NRASOD_1 / D_RASOD_1 / ...
LIM_FITNESS = Cst.L_FIT_MAX         # L_FIT_1 / L_FIT_MAX / L_FIT_UNLIM
NUM_REPETITIONS = 10             # 1 / 10 / 100 / 1000 / 10000
NUM_ITERATIONS = 10000000           # 1 (C_SOLVE_LE) / 10000000 (C_TSTEPPING)
MAX_COUNTER = 5*NUM_REPETITIONS     # 5*NUM_REPETITIONS

POP_SIZE = 0                        # 0 corresponds to an infinite population
NUM_ALL_INI = 10
NUM_PATHO_INI = 200
LOST_THRESHOLD = 1.0E-05
DISREG_THRESHOLD = 1.0E-08
NUM_ALL_THRESHOLD = 100
LAMBDA_IN = 1.0
DELTA_IN = 0.065225
ALPHA_IN = 0.4
UNICL_STEP = 2
NUM_ENV_FIT_CL = 0 # default = 0 -> no environm. fitness distr.
MIN_ALL_FIT = 0.0
MAX_ALL_FIT = 0.1
MIN_GT_FIT = 0.0
MAX_GT_FIT = 1.0
MOD_UPDT = round(NUM_ITERATIONS/100000)
MOD_DISP = round(NUM_ITERATIONS/1000)
MOD_STEP = round(NUM_ITERATIONS/100)
MOD_RATE = round(NUM_ITERATIONS/100)
if CALC_MODE == Cst.C_SOLVE_LE:
    NUM_ITERATIONS = 1
    LOST_THRESHOLD = 10**(-Cst.RPR_15)
    if POP_SIZE > 0:
        LOST_THRESHOLD = 1/float(2*POP_SIZE)
    MOD_UPDT = 1
    MOD_DISP = 1
    MOD_STEP = 1
    MOD_RATE = 1
MOD_FILE = NUM_ITERATIONS
D_FLAGS = {'doCalcs': doCalcs, 'extInf': extractInfo, 'plotMC': plotModelCmp,
           'breakEarly': True, 'delNonHuge': False, 'saveAll': True}
D_CLASS_OUTC = {'wee': -4, 'small': -2, 'normal': 0, 'large': 2, 'huge': 4}
D_N_MODEL = {Cst.D_SOD_1: Cst.N_SOD_S,
             Cst.D_SASOD_1: Cst.N_SAsOD_S,
             Cst.D_RASOD_1: Cst.N_RAsOD_S,
             Cst.D_NRASOD_1: Cst.N_DAA_S}
ND_RESULTS = '_Results'
ND_PLOTS = '_Plots'
NF_STEPRES = 'StepResult'
NF_RATERES = 'RateResult'
NF_GPRES = 'GenepoolResult'
NF_SGPRES = 'ShortSortedGenepoolResult'
NF_HA = 'HeterozygoteAdvantageResult'
NF_HAC = 'HetAdvCorrect'
NF_REPRES = '_RepetitionResult'
NF_BSR = 'BasicResult'
NF_NAS = 'DistrNumAllelesResult'
NF_DFS = 'DistrDiffFitResult'
NF_TKL = 'A_TaskList'
NF_STS = 'B_StepSizes'
NF_PRB = 'C_Problems'
NF_DIC = 'D_Dictionary'
NF_SUM = 'E_Summary'
NF_MCP = 'F_ModelComparison'
NF_MDC = 'G_ModelDistrComp'
NF_NAT = 'NumberOfAlleles'
NF_DFT = 'IntrMeritDifference'
LNF_HACOL = ['ProportionHomozygotes', 'ProportionHeterozygotes',
             'AvFitHomozygotes', 'AvFitHeterozygotes', 'HeterozygoteAdvAbs',
             'HeterozygoteAdvRel']
NFE_TXT = '.dat'
NFE_PLT = '.pdf'

# initiate variables
calcMode = CALC_MODE
inAllFDist = IN_ALL_FDIST
inAllProp = IN_ALL_PROP
dominType = DOM_TYPE
limFitness = LIM_FITNESS
lambdaIn = LAMBDA_IN
deltaIn = DELTA_IN
alphaIn = ALPHA_IN
uniClStep = UNICL_STEP
numEFitCl = NUM_ENV_FIT_CL
minFitn = MIN_ALL_FIT
maxFitn = MAX_ALL_FIT
minGTFitn = MIN_GT_FIT
maxGTFitn = MAX_GT_FIT
lostThr = LOST_THRESHOLD
disregThr = DISREG_THRESHOLD
numAllThr = NUM_ALL_THRESHOLD
numRepetitions = NUM_REPETITIONS
numIters = NUM_ITERATIONS
maxCounter = max(MAX_COUNTER, NUM_REPETITIONS)
popSize = POP_SIZE
numAllI = NUM_ALL_INI
numPatI = NUM_PATHO_INI
modUpdt = MOD_UPDT
modDisp = MOD_DISP
modStep = MOD_STEP
modRate = MOD_RATE
modFile = MOD_FILE
dFlags = D_FLAGS
dClOutc = D_CLASS_OUTC
dNModel = D_N_MODEL
lModels = [Cst.D_NRASOD_1, addModel]
lNModel = [dNModel[cDT] for cDT in lModels]

# variables for timings
startTimeSim = 0.0
endTimeSim = 0.0

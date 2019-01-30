################################################################################
# GeneEvolutionMainV13.py #
################################################################################

import random, time
from Input import *
from Constants import C_SOLVE_LE, C_TSTEPPING, RPR_04
from SpecialFunctions import (genDictInp, buildTaskLD, genDictRepRes, incCTask,
                              getFinResDicts, getLModels, updateDTPR,
                              genDictAllFit, genDictAllProp,
                              genDictGenotypeFit, genGenepoolDict, 
                              oneTimeStepDynamics, updateDicts, calcFinalDicts)
from GenericFunctions import checkIfSuccess
from WriteAndPrint import (writeLTasks, deleteAllFiles, printCurTask,
                           printCurRepet, printTaskFinished, openOutFiles,
                           writeResults, writeDistsRes, writeBasicRes,
                           writeHetAdvCorrectRes, writeRepRes,
                           printElapsedTimeSim)
from GeneEvolutionTools import extractAllInfo

print('GeneEvolutionV13')
print('Determines the distribution of allele frequencies at equilibrium.')

# MAIN PROGRAM #
startTimeSim = time.time()
print('+'*20 + ' START', time.ctime(startTimeSim), '+'*20)
cWD = os.getcwd()
print('Current working directory:', cWD)
dictInp = genDictInp(calcMode, inAllFDist, inAllProp, limFitness, dominType,
                     lambdaIn, deltaIn, alphaIn, uniClStep, numEFitCl, minFitn,
                     maxFitn, minGTFitn, maxGTFitn, lostThr, disregThr,
                     numAllThr, numRepetitions, numIters, maxCounter, popSize,
                     numAllI, numPatI, modUpdt, modDisp, modStep, modRate,
                     modFile, dFlags, dClOutc, dNModel, lModels, lNModel,
                     ND_RESULTS, ND_PLOTS, NF_STEPRES, NF_RATERES, NF_GPRES,
                     NF_SGPRES, NF_HA, NF_HAC, NF_REPRES, NF_BSR, NF_NAS,
                     NF_DFS, NF_TKL, NF_STS, NF_PRB, NF_DIC, NF_SUM, NF_MCP,
                     NF_MDC, NF_NAT, NF_DFT, LNF_HACOL, NFE_TXT, NFE_PLT)
lDTasks = buildTaskLD(lDTasksIn, taskMode)
nModels, nTasks = len(dictInp['lModels']), len(lDTasks)
writeLTasks(dictInp, lDTasks, cWD)
random.seed()
if dictInp['dFlags']['doCalcs']:
    for iTask, cDTask in enumerate(lDTasks):
        incCTask(dictInp, cDTask)
        printCurTask(iTask, nTasks, dictInp['cDTask'])
        dictRepRes = genDictRepRes(dictInp)
        dictDistNAll, dictDistDFit, dictBasicRes = getFinResDicts(dictInp)
        curRep, curCt = 1, 1
        deleteAllFiles(dictInp, cWD)
        while curRep <= dictInp['numRep'] and curCt <= dictInp['maxCounter']:
            print('Current folder:', cWD.split('/')[-1])
            printCurRepet(dictInp, 0, iTask, nModels, nTasks, curRep)
            dictAFit = genDictAllFit(dictInp)
            lFCnv, lFRnk = [False]*nModels, [False]*nModels
            lCModels = getLModels(dictInp, curRep)
            for iModel, cModel in enumerate(lCModels):
                cNModel = updateDTPR(dictInp, cModel)
                printCurTask(iTask, nTasks, dictInp['cDTask'], cNModel)
                dictAProp = genDictAllProp(dictInp, dictAFit)
                dictGenType = genDictGenotypeFit(dictInp, dictAFit)
                dictGeneP, dictRes = genGenepoolDict(dictInp, dictAFit,
                                                     dictAProp, dictGenType)
                if dictInp['dFlags']['saveAll']:
                    openOutFiles(dictInp, dictGeneP, curRep, cWD)
                isC, isFRnk = False, False
                if dictInp['calcMd'] == C_SOLVE_LE:
                    goOn = writeResults(dictInp, dictGeneP, dictGenType,
                                        dictRes, dictRepRes, curRep, 0, cWD,
                                        isC)
                    isC, lFRnk[iModel] = calcFinalDicts(dictInp, dictGeneP,
                                                        dictGenType, dictRes)
                    if lFRnk[iModel]:
                        goOn = writeResults(dictInp, dictGeneP, dictGenType,
                                            dictRes, dictRepRes, curRep,
                                            numIters, cWD, isC)
                elif dictInp['calcMd'] == C_TSTEPPING:
                    for timeStep in range(numIters + 1):
                        if timeStep > 0 and not isC:
                            isC = oneTimeStepDynamics(dictInp, dictGeneP,
                                                      dictGenType, dictRes)
                        goOn = writeResults(dictInp, dictGeneP, dictGenType,
                                            dictRes, dictRepRes, curRep,
                                            timeStep, cWD, isC)
                        if not goOn:
                            break
                else:
                    print('ERROR: Calculation mode not correctly specified.')
                updateDicts(dictInp, dictDistNAll, dictDistDFit, dictBasicRes,
                            dictGeneP, dictRes, dictRepRes, lFRnk[iModel],
                            curRep)
            curCt += 1
            if (checkIfSuccess(lFRnk, len(lCModels)) or
                dictInp['calcMd'] == C_TSTEPPING):
                curRep += 1
                printElapsedTimeSim(startTimeSim, time.time())
#         writeDistsRes(dictInp, dictDistNAll, dictDistDFit, cWD)
        writeBasicRes(dictInp, dictBasicRes, cWD)
        writeHetAdvCorrectRes(dictInp, dictRepRes, cWD)
        writeRepRes(dictInp, dictRepRes, cWD)
        printTaskFinished(iTask, nTasks)
if dictInp['dFlags']['extInf']:
    extractAllInfo(dictInp, cWD, RPR_04)
printElapsedTimeSim(startTimeSim, time.time(), 'Total time')
print('*'*20 + ' DONE', time.ctime(time.time()), '*'*20)

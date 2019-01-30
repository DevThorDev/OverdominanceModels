################################################################################
# SpecialFunctions.py #
################################################################################

import random, math

import numpy as np
from scipy import linalg

import Constants as Cst
from GenericFunctions import (addToCountDict, convSubSetToBinChain,
                              randPartitionRange, calcRelFitness, calcMean,
                              getKeyOfMaxDictVal, getNElD, addNewD, modLIdx,
                              safeRemoveEl, checkIfSqMat, removeRowColMat)

# ### FUNCTIONS ################################################################
def getDictNumRep(dNMdl, lNMdl, allFDst, nRepStd = 1, nRepSpec = 1):
    dNRep = {cNMdl: nRepStd for cNMdl in lNMdl}
    if allFDst in [Cst.A_FDIST_UNIFORM]:
        for cNMdl in lNMdl:
            if cNMdl not in [dNMdl[Cst.D_NRASOD_1]]:
                dNRep[cNMdl] = nRepSpec
    return dNRep

def genDictInp(calcMd, allFDist, allProp, limFit, domType, lambdaIn, deltaIn,
               alphaIn, uniClStep, numEFitCl, minAFit, maxAFit, minGTFit,
               maxGTFit, lostThr, disregThr, numAllThr, numRep, numIter, maxCt,
               popSz, numAllIni, numPathoIni, modUpdt, modDisp, modStep,
               modRate, modFile, dFlags, dClOutc, dNModel, lModels, lNModel,
               nD_Results, nD_Plots, nF_StepRes, nF_RateRes, nF_GPRes,
               nF_SGPRes, nF_HetAdv, nF_HAC, nF_RepRes, nF_BSR, nF_NAS, nF_DFS,
               nF_TKL, nF_STS, nF_PRB, nF_Dic, nF_Sum, nF_MCP, nF_MDC, nF_NAT,
               nF_DFT, lNF_HACol, nFE_Txt, nFE_Plt):
    theDict = {'calcMd': calcMd,
               'allFDist': allFDist,
               'allProp': allProp,
               'limFit': limFit,
               'domType': domType,
               'lambdaIn': lambdaIn,
               'deltaIn': deltaIn,
               'alphaIn': alphaIn,
               'uniClStep': uniClStep,
               'numEFitCl': numEFitCl,
               'minAFit': minAFit,
               'maxAFit': maxAFit,
               'minGTFit': minGTFit,
               'maxGTFit': maxGTFit,
               'lostThr': lostThr,
               'disregThr': min(disregThr, lostThr),
               'numAllThr': numAllThr,
               'numRep': numRep,
               'numIter': numIter,
               'maxCounter': maxCt,
               'popSz': popSz,
               'numAllIni': numAllIni,
               'numPathoIni': numPathoIni,
               'modUpdt': modUpdt,
               'modDisp': modDisp,
               'modStep': modStep,
               'modRate': modRate,
               'modFile': modFile,
               'dNumRep': getDictNumRep(dNModel, lNModel, allFDist, numRep),
               'dFlags': dFlags,
               'dClOutc': dClOutc,
               'dNModel': dNModel,
               'lModels': lModels,
               'lNModel': lNModel,
               'nD_ResultsC': nD_Results,
               'nD_Results': nD_Results,
               'nD_Plots': nD_Plots,
               'nF_StepRes': nF_StepRes,
               'nF_RateRes': nF_RateRes,
               'nF_GPRes': nF_GPRes,
               'nF_SGPRes': nF_SGPRes,
               'nF_HetAdv': nF_HetAdv,
               'nF_HAC': nF_HAC,
               'nF_RepRes': nF_RepRes,
               'nF_BSR': nF_BSR,
               'nF_NAS': nF_NAS,
               'nF_DFS': nF_DFS,
               'nF_TKL': nF_TKL,
               'nF_STS': nF_STS,
               'nF_Prb': nF_PRB,
               'nF_Dic': nF_Dic,
               'nF_Sum': nF_Sum,
               'nF_MCP': nF_MCP,
               'nF_MDC': nF_MDC,
               'nF_NAT': nF_NAT,
               'nF_DFT': nF_DFT,
               'lNF_HACol': lNF_HACol,
               'nFE_Txt': nFE_Txt,
               'nFE_Plt': nFE_Plt}
    return theDict

def buildTaskLD(lDIn, tMode = 'Mult'):
    lDOut = lDIn
    if tMode == 'Mult':     # this is currently the only mode implemented
        lDOut = []
        for dIn in lDIn:
            nElTot, lNElCat, nCat = getNElD(dIn)
            lIs, lDOutC = [0]*nCat, [{} for _ in range(nElTot)]
            iElTot = 0
            while lIs != [lNElCat[k] - 1 for k in range(len(lNElCat))]:
                addNewD(dIn, lDOutC, lIs, nCat, iElTot)
                modLIdx(lIs, lNElCat, nCat)
                iElTot += 1
            addNewD(dIn, lDOutC, lIs, nCat, iElTot)
            lDOut += lDOutC
    print('+'*32, 'Number of tasks:', len(lDOut), '+'*32)
    return lDOut

def genDictRepRes(dictI):
    dictRpRs = {cNMdl: {} for cNMdl in dictI['lNModel']}
    for cNMdl in dictI['lNModel']:
        dClOutcInc = {}
        for oneCl in dictI['dClOutc']:
            dClOutcInc[oneCl] = 0
        dictRpRs[cNMdl] = {'Problems': {'NotConv': [], 'SumHetHom': {},
                                        'SumAProp': {}, 'FSingular': []},
                           'clOutcInc': dClOutcInc,
                           'dAllPers': {},
                           'dMaxNA': {},
                           'dHetAdvCorrect': {}}
    return dictRpRs

def incCTask(dInp, cDT):
    if cDT['maxAFit'] >= cDT['minAFit']:
        dInp['minAFit'] = cDT['minAFit']
        dInp['maxAFit'] = cDT['maxAFit']
    else:
        dInp['minAFit'] = cDT['maxAFit']
        dInp['maxAFit'] = cDT['minAFit']
    cDT['maxAFit'] = dInp['maxAFit']
    cDT['minAFit'] = dInp['minAFit']
    rngAFit = dInp['maxAFit'] - dInp['minAFit']
    assert rngAFit >= 0
    dInp['numAllIni'] = cDT['numAllIni']
    dInp['numPathoIni'] = cDT['numAllIni']*cDT['multNumPatho']
    if rngAFit > 0:
        dInp['numPathoIni'] = round(dInp['numPathoIni']/rngAFit)
    dInp['uniClStep'] = cDT['multNumPatho']
    dInp['cDTask'] = cDT

def getFinResDicts(dictI):
    dictDstNA, dictDstDF, dictBasRs = {}, {}, {}
    for cNMdl in dictI['lNModel']:
        dictDstNA[cNMdl] = {}
        dictDstDF[cNMdl] = {}
        dictBasRs[cNMdl] = {}
    return dictDstNA, dictDstDF, dictBasRs

def getLModels(dictI, cRep):
    lCMdl = dictI['lModels']
    if dictI['allFDist'] in [Cst.A_FDIST_UNIFORM]:
        if cRep > 1:
            if Cst.D_NRASOD_1 in dictI['lModels']:
                lCMdl = [Cst.D_NRASOD_1]
            else:
                lCMdl = []
    return lCMdl

def getPRes(dInp, domTp = -1):
    if domTp == -1:     # not specified - take domTp from input dictionary
        domTp = dInp['domType']
    pResExt = '_'  + dInp['calcMd'] + '_' + dInp['dNModel'][domTp]
    for cK, cV in sorted(dInp['cDTask'].items()):
        pResExt += ('__' + str(cK) + '_' + str(cV).replace('.', 'p'))
    return pResExt

def updateDTPR(dInp, cMdl):
    dInp['domType'] = cMdl
    dInp['nD_Results'] = dInp['nD_ResultsC'] + getPRes(dInp)
    return dInp['dNModel'][cMdl]

def genDictFDistUniformCl(dI):
    nEl = dI['numPathoIni']
    cStart, cStep = int(dI['uniClStep']/2), int(dI['uniClStep'])
    lAll = []
    while len(lAll) < dI['numAllIni']:
        for n1s in range(math.ceil(nEl*dI['minAFit']) + cStart,
                         math.floor(nEl*dI['maxAFit']) + cStart, cStep):
            lCur = random.sample(range(nEl), n1s)
            lAll.append(lCur)
            if len(lAll) >= dI['numAllIni']:
                break
    return convSubSetToBinChain(lAll, nEl)

def genDictFDistRandom(dI):
    nEl = dI['numPathoIni']
    dAll = {}
    for i in range(dI['numAllIni']):
        offsAFit = dI['uniClStep']/(2*nEl)
        minAFit, maxAFit = dI['minAFit'] + offsAFit, dI['maxAFit'] - offsAFit
        n1s = round(nEl*np.random.uniform(minAFit, maxAFit))
        lRespPat = [0]*(nEl - n1s) + [1]*n1s
        lRespPat = np.random.permutation(lRespPat)
        dAll[i] = [lRespPat.tolist(), n1s/nEl]
    return dAll

def genDictAllFit(dictI):
    dAF = {}
    if dictI['allFDist'] == Cst.A_FDIST_UNIFORM:    # "uniform" [minAFit, maxAFit]
        dAF = genDictFDistUniformCl(dictI)
    elif dictI['allFDist'] == Cst.A_FDIST_RANDOM:   # random from uniform distribution
        dAF = genDictFDistRandom(dictI)
    else: pass
    return dAF

def genDictAllProp(dictI, dictAF):
    numAll, dictAP = len(dictAF), {}
    if dictI['popSz'] > 0:
        numAPop = 2*dictI['popSz']
        if dictI['allProp'] == Cst.A_PROP_RANDOM:
            lIAllDr = np.random.randint(numAll, size = numAPop)
            for j in range(numAPop):
                addToCountDict(dictAP, lIAllDr[j], 1.0/numAPop)
        elif dictI['allProp'] == Cst.A_PROP_UNIFORM:
            for j in range(numAll):
                dictAP[j] = float((numAPop//numAll)/numAPop)
            lIAllDr = np.random.randint(numAll, size = numAPop%numAll)
            for j in range(numAPop%numAll):
                addToCountDict(dictAP, lIAllDr[j], 1.0/numAPop)
        else: pass
    else:
        if dictI['allProp'] == Cst.A_PROP_RANDOM:
            lProps = randPartitionRange(numAll)
            for j in range(numAll):
                dictAP[j] = lProps[j]
        elif dictI['allProp'] == Cst.A_PROP_UNIFORM:
            for j in range(numAll):
                dictAP[j] = round(1.0/numAll, Cst.RPR_15)
        else: pass
    return dictAP

def getDefGTDict(dI, dAF):
    dGTF = {}
    for aC1 in dAF:
        for aC2 in dAF:
            if aC2 <= aC1:
                [fPGT, ovGT] = calcRelFitness(dAF[aC1][0], dAF[aC2][0])
                fitGT = dI['minGTFit'] + (dI['maxGTFit'] - dI['minGTFit'])*fPGT
                dGTF[(aC1, aC2)] = [fitGT, ovGT]
    return dGTF

def G2FitNoHA(dI, fC1, fC2):
    fGT = dI['minGTFit']
    ovGT = max(fC1, fC2) - min(fC1, fC2)
    if dI['domType'] == Cst.D_RND_DOM:
        # draw either first or second chromosome of climate fitness gene
        if random.randint(1, 2) == 1:
            fGT = fC1
        else:
            fGT = fC2
    elif dI['domType'] == Cst.D_BETTER_AD_DOM:
        fGT = max(fC1, fC2)
    elif dI['domType'] == Cst.D_WORSE_AD_DOM:
        fGT = min(fC1, fC2)
    elif dI['domType'] == Cst.D_CODOM:
        fGT = round(calcMean([fC1, fC2]), Cst.RPR_15)
    else: pass
    return [fGT, ovGT]

def G2FitWithHA(dI, aC1, aC2, fC1, fC2):
    fGT = dI['minGTFit']
    ovGT = max(fC1, fC2) - min(fC1, fC2)
    if dI['domType'] == Cst.D_RASOD_1:
        if aC1 == aC2:      # homozygote
            fGT = fC1
        else:               # heterozygote
            fGT = round(fC1 + (1 - dI['lambdaIn']*fC1)*fC2, Cst.RPR_15)
        ovGT = dI['lambdaIn']*fC1*fC2
    elif dI['domType'] == Cst.D_SASOD_0_MX:
        if aC1 == aC2:      # homozygote
            fGT = fC1
        else:               # heterozygote
            fGT = round(max(fC1, fC2) + dI['deltaIn'], Cst.RPR_15)
    elif dI['domType'] == Cst.D_SASOD_1:
        if aC1 == aC2:      # homozygote
            fGT = fC1
        else:               # heterozygote
            fGT = round(calcMean([fC1, fC2]) + dI['deltaIn'], Cst.RPR_15)
    elif dI['domType'] == Cst.D_SOD_1:
        if aC1 == aC2:      # homozygote
            fGT = fC1
        else:               # heterozygote
            fGT = dI['maxGTFit']
    else: pass
    return [fGT, ovGT]

def convGene2Fit(dI, dAF, allC1, allC2):
    fC1 = dAF[allC1][1]
    fC2 = dAF[allC2][1]
    if dI['domType'] in [Cst.D_RND_DOM, Cst.D_BETTER_AD_DOM,
                         Cst.D_WORSE_AD_DOM, Cst.D_CODOM]:
        [fGT, ovGT] = G2FitNoHA(dI, fC1, fC2)
    else:
        [fGT, ovGT] = G2FitWithHA(dI, allC1, allC2, fC1, fC2)
    # cap genotype fitness and add excess to overlap
    fGT = max(fGT, dI['minGTFit'])
    if dI['limFit'] == Cst.L_FIT_1 or dI['limFit'] == Cst.L_FIT_MAX:
        if fGT > dI['maxGTFit']:
            ovGT += (fGT - dI['maxGTFit'])
        fGT = min(fGT, dI['maxGTFit'])
    return [fGT, ovGT]

def getGenericGTDict(dI, dAF):
    dGTF = {}
    for aC1 in dAF:
        for aC2 in dAF:
            if aC2 <= aC1:
                dGTF[(aC1, aC2)] = convGene2Fit(dI, dAF, aC1, aC2)
    return dGTF

def genDictGenotypeFit(dictI, dictAF):
    dictGTF = {}
    if dictI['domType'] in [Cst.D_NRASOD_1]:
        dictGTF = getDefGTDict(dictI, dictAF)
    else:
        dictGTF = getGenericGTDict(dictI, dictAF)
    return dictGTF

def removeLowFreqAlleles(dI, dGP, dR, lIAll):
    for iAll in lIAll:
        if dGP[iAll][2] < dI['lostThr'] and iAll in dR['lAllPers']:
            dR['lAllPers'].remove(iAll)
            dR['nAllTSLost'] += 1
            dR['nAllTotLost'] += 1
            safeRemoveEl(dGP, iAll, 2, dGP[iAll][2])

def updateResDict(dI, dGP, dR, popF = -1.0):
    dR['propZyg']['propHomoZyg'] = 0.0
    nAllDisreg = len(dR['lAllPers']) - len(dGP)
    dR['nAllTSLost'] += nAllDisreg
    dR['nAllTotLost'] += nAllDisreg
    dR['lAllPers'], lIAll = list(dGP), list(dGP)
    removeLowFreqAlleles(dI, dGP, dR, lIAll)
    lIAll = list(dGP)
    for curIAll in lIAll:
        dR['propZyg']['propHomoZyg'] += dGP[curIAll][2]*dGP[curIAll][2]
    dR['nAllPers'] = len(dR['lAllPers'])
    dR['propZyg']['propHeteroZyg'] = 1.0 - dR['propZyg']['propHomoZyg']
    if popF >= 0:   # only set if a new population fitness value is given 
        dR['popFitness'] = round(popF, Cst.RPR_15)

def calcFitInc_FDst(curA, othA, dI, dAP, dGTF):
    mFInc = 0
    if dI['numEFitCl'] > 0:
        pass
    else:
        mFInc = dAP[othA]*dGTF[(max(curA, othA), min(curA, othA))][0]
    return mFInc

def calcMargFit(curA, dI, dAP, dGTF, roundPrec = 16):
    mF = 0.0
    for othA in dAP:
        mF += calcFitInc_FDst(curA, othA, dI, dAP, dGTF)
    return round(mF, roundPrec)

def genGenepoolDict(dictI, dictAF, dictAP, dictGTF):
    popFit = 0.0
    dictGeneP = {}
    dictZyg = {'propHomoZyg': 0.0, 'propHeteroZyg': 1.0}
    dictR = {'lAllPers': list(dictAF), 'nAllPers': len(dictAF),
             'nAllTSLost': 0, 'nAllTotLost': 0, 'propZyg': dictZyg}
    for curIAll in dictAF:
        margF = calcMargFit(curIAll, dictI, dictAP, dictGTF, Cst.RPR_15)
        if curIAll in dictAP:
            dictGeneP[curIAll] = [dictAF[curIAll][0], dictAF[curIAll][1],
                                  dictAP[curIAll], margF]
            popFit += dictAP[curIAll]*margF
        else:
            print('Error: keys of dictionaries do not match')
            dictGeneP[curIAll] = [0.0]*4
    updateResDict(dictI, dictGeneP, dictR, popFit)
    return dictGeneP, dictR

def updateDictGeneP(dictI, dictAP, dictGP, dictGTF, dictR):
    popFitNew = 0.0
    minMargFit, maxMargFit = dictI['maxGTFit'], dictI['minGTFit']
    for cIA in list(dictGP):
        if cIA in dictAP:
            newMargFit = calcMargFit(cIA, dictI, dictAP, dictGTF, Cst.RPR_15)
            dictGP[cIA] = dictGP[cIA][:2] + [dictAP[cIA], newMargFit]
            popFitNew += dictAP[cIA]*newMargFit
            if newMargFit < minMargFit:
                minMargFit = newMargFit
            if newMargFit > maxMargFit:
                maxMargFit = newMargFit
        else:
            del dictGP[cIA]
    # population fitness is the only result value needed - hence recalculated
    dictR['popFitness'] = round(popFitNew, Cst.RPR_15)
    return minMargFit, maxMargFit

def oneTimeStepDynamics(dictI, dictGeneP, dictGTF, dictR):
    dictAP, lIA = {}, list(dictGeneP)
    # create dictAP
    for cIA in lIA:
        curProp = dictGeneP[cIA][2]*dictGeneP[cIA][3]/dictR['popFitness']
        if curProp >= dictI['disregThr']:
            dictAP[cIA] = curProp
        else:
            safeRemoveEl(dictGeneP, cIA, 2, dictGeneP[cIA][2])
    tMargFit = updateDictGeneP(dictI, dictAP, dictGeneP, dictGTF, dictR)
    return (tMargFit[1] - tMargFit[0])/tMargFit[0]*100.0 < 10**(-Cst.RPR_11)

def getMostCommonAllele(dGP, lAllPs, lAllUs):
    dictIdx = 2
    outString = 'Error: No allele has a positive frequency...'
    return getKeyOfMaxDictVal(dGP, lAllPs, lAllUs, dictIdx, outString)

def getFittestAllele(dGP, lAllPs, lAllUs):
    dictIdx = 1
    outString = 'Error: No allele has a positive fitness...'
    return getKeyOfMaxDictVal(dGP, lAllPs, lAllUs, dictIdx, outString)

def updateDicts(dI, dDstNA, dDstDF, dBasR, dGP, dR, dRR, isFR, cRp):
    fitMin, fitMax, sumProp = 1.0, 0.0, 0.0
    for cAll in dR['lAllPers']:
        cFit, cProp = dGP[cAll][1], dGP[cAll][2]
        if cFit < fitMin:
            fitMin = cFit
        if cFit > fitMax:
            fitMax = cFit
        sumProp += cProp
    diffFitTot = round(fitMax - fitMin, Cst.RPR_14)
    cNMod = dI['dNModel'][dI['domType']]
    if dR['nAllPers'] > 0 and diffFitTot >= 0:
        addToCountDict(dDstNA[cNMod], dR['nAllPers'])
        addToCountDict(dDstDF[cNMod], diffFitTot)
        dBasR[cNMod][cRp] = [dR['nAllPers'], diffFitTot]
        if isFR and round(sumProp, Cst.RPR_12) != 1.0:
            dRR[cNMod]['Problems']['SumAProp'][cRp] = round(sumProp, Cst.RPR_12)

def convertDictGTF(dGTF):
    nAll = max([max(cK) for cK in dGTF]) + 1
    F = np.zeros((nAll, nAll))
    for cIAC1 in range(nAll):
        for cIAC2 in range(nAll):
            if (cIAC1, cIAC2) in dGTF:
                cGTF = dGTF[(cIAC1, cIAC2)][0]
                F[cIAC1, cIAC2] = cGTF
                F[cIAC2, cIAC1] = cGTF
    return F

def checkMatFullRank(sqMat):
    hasFullRank = True
    if np.linalg.matrix_rank(sqMat) < sqMat.shape[0]:
        print('WARNING: Genotype fitness matrix does not have full rank.')
        hasFullRank = False
    return hasFullRank

def adaptIndexMap(lMapIA, lIRet):
    for iI, iRt in enumerate(lIRet):
        lMapIA[iI] = lMapIA[iRt]
    del lMapIA[len(lIRet):]

def getVAllProp(F):
    vAP = np.zeros(F.shape[0])
    if checkMatFullRank(F):
        isFR = True
        vAP = np.linalg.solve(F, np.ones(F.shape[0]))
        vAP /= sum(vAP)
    else:
        isFR = False
    return isFR, vAP
 
def reduceFitMat(F, lstT, lMapIA):
    isFR, vAP = getVAllProp(F)
    if isFR:
        lIRet = []
        for iAP, cAP in enumerate(vAP):
            if cAP >= lstT:
                lIRet.append(iAP)
        adaptIndexMap(lMapIA, lIRet)
        F = removeRowColMat(F, lIRet)
        isFR, vAP = getVAllProp(F)
    return F, isFR, vAP

def calcAllProp(F, lstThr):
    dAP, isFRk = {}, False
    if checkIfSqMat(F):
        vAP, lMapIA = np.zeros(F.shape[0]), list(range(F.shape[0]))
        while min(vAP) < lstThr:
            F, isFRk, vAP = reduceFitMat(F, lstThr, lMapIA)
            if not isFRk:
                break
        for iV, iA in enumerate(lMapIA):
            dAP[iA] = vAP[iV]
    return dAP, isFRk

def calcFinalDicts(dictI, dictGeneP, dictGTF, dictR):
    lostThr = max(dictI['lostThr'], 10**(-Cst.RPR_15))
    dictAP, isFRk = calcAllProp(convertDictGTF(dictGTF), lostThr)
    tMargFit = updateDictGeneP(dictI, dictAP, dictGeneP, dictGTF, dictR)
    percDevMargFit = 0.0
    if tMargFit[0] > 0:
        percDevMargFit = (tMargFit[1] - tMargFit[0])/tMargFit[0]*100.0
        if percDevMargFit >= 10**(-Cst.RPR_11):
            print('WARNING: Marginal fitness deviation:', percDevMargFit, '%')
    else:
        print('WARNING: Minimal marginal fitness is', tMargFit[0])
    return (percDevMargFit < 10**(-Cst.RPR_11)), isFRk

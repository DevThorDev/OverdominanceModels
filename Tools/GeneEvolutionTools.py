import os, math
import numpy as np
import matplotlib.pyplot as plt

from Constants import RPR_12
from GenericFunctions import createDir, addToCountDict

# ### FUNCTIONS ### ------------------------------------------------------------
def getKeyInfo(nD):
    lK = ['', ['']*4]
    lSplC = nD.split('__')
    if len(lSplC) >= 1:
        lSpl0 = lSplC[0].split('_')
        lSpl1 = lSplC[1:]
        if len(lSpl0) >= 4:
            lK[0] = lSpl0[3]
        if len(lSpl1) >= 4:
            for cI, cStr in enumerate(lSpl1):
                lSpl2 = cStr.split('_')
                if len(lSpl2) >= 2:
                    lK[1][cI] = float(lSpl2[1].replace('p', '.'))
                    if cI >= 2:
                        lK[1][cI] = int(lK[1][cI])
            lK[1] = tuple(lK[1])
    return tuple(lK)

def addValToDictInf(dI, cK, tK, cV, isList = False):
    if tK[0] not in dI[cK]:
        dI[cK][tK[0]] = {}
    if isList:
        if tK[1] not in dI[cK][tK[0]]:
            dI[cK][tK[0]][tK[1]] = [cV]
        else:
            dI[cK][tK[0]][tK[1]].append(cV)
    else:
        dI[cK][tK[0]][tK[1]] = cV

def getVal(dI, lSp, cK, tK, cvDtInt = True):
    cPos = 0
    if cK == 'lFreqs':
        cPos = 1
        cvDtInt = True
    cV = float(lSp[cPos])
    if cvDtInt:
        cV = int(round(cV))
    addValToDictInf(dI, cK, tK, cV, True)
    return cV

def calcMeanNoNA(lVs):
    if len(lVs) > 0:
        cM, nVs = 0.0, 0
        for cV in lVs:
            if cV != 'NA':
                nVs += 1
                cM += cV
        if nVs > 0:
            return cM/nVs
        else:
            return 'NA'
    else:
        return 'NA'

def calcWMean(lVs, lFrs):
    cWM = 0.0
    if len(lVs) == len(lFrs):
        for k in range(len(lVs)):
            cWM += lVs[k]*lFrs[k]
    else:
        print('length of values =', len(lVs),
              'but unequal to length of frequencies =', len(lFrs))
    return cWM/sum(lFrs)

def extractDistributions(dInp, tD, cD, nFBR, nFNA, nFDF):
    dDstNA, dDstDF = {}, {}
    fBR = open(os.path.join(tD, cD, nFBR), 'r')
    for cI, cL in enumerate(fBR):
        if cI > 0:
            cRp, cNA, cDF = cL[:-1].split(' ')
            cRp, cNA, cDF = int(cRp), int(cNA), float(cDF)
#             print('cRp:', cRp, '- cNA:', cNA, '- cDF:', cDF)
            addToCountDict(dDstNA, cNA)
            addToCountDict(dDstDF, round(cDF, RPR_12))
    fBR.close()
    # open number of alleles distribution file
    fRes = open(os.path.join(tD, cD, nFNA), 'w')
    fRes.write('NumberAlleles TimesRecorded\n')
    for cKey in sorted(dDstNA, reverse = True):
        fRes.write(str(cKey) + ' ' + str(dDstNA[cKey]) + '\n')
    fRes.close()
    # open fitness difference distribution file
    fRes = open(os.path.join(tD, cD, nFDF), 'w')
    fRes.write('FitnessDifferenceAbs TimesRecorded\n')
    for cKey in sorted(dDstDF, reverse = True):
        fRes.write(str(cKey) + ' ' + str(dDstDF[cKey]) + '\n')
    fRes.close()

def extractValsFreqs(dInf, tD, cD, nF, cnvDtInt = True):
    tK = getKeyInfo(cD)
    pF = os.path.join(tD, cD, nF)
    if os.path.isfile(pF):
        lVals, lFreqs = [], []
        fS = open(pF, 'r')
        for cI, cL in enumerate(fS):
            if cI >= 1 and len(cL) >= 1:
                lSpl = cL.split(' ')
                if len(lSpl) >= 2:
                    lVals.append(getVal(dInf, lSpl, 'lVals', tK, cnvDtInt))
                    lFreqs.append(getVal(dInf, lSpl, 'lFreqs', tK, cnvDtInt))
        fS.close()
        assert len(lVals) > 0 and len(lFreqs) > 0 and len(lVals) == len(lFreqs)
        addValToDictInf(dInf, 'wMean', tK, calcWMean(lVals, lFreqs))
        addValToDictInf(dInf, 'maxV', tK, max(lVals))
        addValToDictInf(dInf, 'minV', tK, min(lVals))

def fillDArrVals(dArrVs, cI, cL, nC, nCDt, nF):
    doBreak = False
    if cI == 0:
        nC = len(cL.split(' '))
        if nC != nCDt + 1:
            print('Error: File', nF, 'has', nC, 'columns!')
            doBreak = True
    else:
        if len(cL) >= 1:
            lSpl = cL.split(' ')
            lSpl[-1] = lSpl[-1][:-1]
            if len(lSpl) == nC:
                for k in range(nCDt):
                    if lSpl[k + 1] != 'NA':
                        dArrVs[k].append(float(lSpl[k + 1]))
            else:
                print('Error: File', nF, 'has', nC, 'columns.')
    return nC, doBreak

def fillDMeanStdDev(dArrVs, dMnSD, nCDt):
    for k in range(nCDt):
        if len(dArrVs[k]) > 0:
            arrVals = np.array(dArrVs[k], dtype = np.float64)
            dMnSD[k] = (np.mean(arrVals), np.std(arrVals))
        else:
            dMnSD[k] = ('NA', 'NA')

def extractOtherData(dODat, tD, cD, nF, dNRep):
    tK = getKeyInfo(cD)
    pF = os.path.join(tD, cD, nF)
    if os.path.isfile(pF):
        nC, nCDt = 0, 6
        dArrVals, dMnSD = {k: [] for k in range(nCDt)}, {}
        fS = open(pF, 'r')
        for cI, cL in enumerate(fS):
            nC, errStop = fillDArrVals(dArrVals, cI, cL, nC, nCDt, nF)
            if errStop:
                break
        fS.close()
        fillDMeanStdDev(dArrVals, dMnSD, nCDt)
        if tK[0] not in dODat['HetAdvData']:
            dODat['HetAdvData'][tK[0]] = {}
        dODat['HetAdvData'][tK[0]][tK[1]] = dMnSD

def getActionIncCt(cL, cS, iC):
    iA = 0              # iA: indicator of action (default: 0)
    if len(cL) >= len(cS):
        if cL[:len(cS)] == '-'*len(cS):
            iA = 1      # a new section begins
            iC += 1     # shows which section we're in
        elif cL[:len(cS)] == cS:
            iA = 2      # first line of new section
        else:
            iA = 3      # line with information
    else:
        if len(cL) >= 1:
            iA = 3      # line with information
    return iA, iC

def getListInfo(cL):
    lI, lSpl = [], cL.split(',')
    for cI, cStr in enumerate(lSpl):
        if cI == 0:
            if len(cStr) >= 1:
                cStr = cStr[1:]
        if cI == len(lSpl) - 1:
            if len(cStr) >= 2:
                cStr = cStr[:-2]
        lI.append(int(cStr.strip()))
    return sorted(lI)

def addToDictInfo(cL, cD):
    lSpl = cL.split(':')
    if len(lSpl) == 2:
        cD[int(lSpl[0])] = float(lSpl[-1].strip())

def useDictInfAsValue(dI, dSum, keyDSum):
    if len(dI) > 0:     # reading info into dictionary is finished
        dSum[keyDSum] = dI
        dI = {}

def fillSubDicts(dNCv, dSHH, dSAF, dCt2, dCt3, kL3, cL, iA, iC):
    if iA == 3:         # if it is a line with information
        if iC == 1:     # 'not converged'
            dNCv[kL3] = getListInfo(cL)
        elif iC == 2:   # 'sum allele frequencies'
            addToDictInfo(cL, dCt2)
        elif iC == 3:   # 'sum frequencies homo- and heterozygotes'
            addToDictInfo(cL, dCt3)
    elif iA == 1:       # a new section begins
        useDictInfAsValue(dCt2, dSHH, kL3)    # add dCt2 as new value to dSHH
        useDictInfAsValue(dCt3, dSAF, kL3)    # add dCt3 as new value to dSAF
        dCt2, dCt3 = {}, {}

def addToDProb(dInf, dProb, kL1, kL2, kL3):
    if kL2 not in dProb[kL1]:
        dProb[kL1][kL2] = {}
    for kL3, cV in dInf.items():
        dProb[kL1][kL2][kL3] = cV

def extractProblems(dProb, tD, cD, nF):
    tK = getKeyInfo(cD)
    pF = os.path.join(tD, cD, nF)
    dNotCnv, dSumHH, dSumAP = {}, {}, {}
    sRp = 'Repetitions'
    if os.path.isfile(pF):
        ct, dCt2, dCt3 = 0, {}, {}
        fS = open(pF, 'r')
        for cL in fS:
            iA, ct = getActionIncCt(cL, sRp, ct)
            fillSubDicts(dNotCnv, dSumHH, dSumAP, dCt2, dCt3, tK[1],
                         cL, iA, ct)
            if ct >= 4 and iA >= 3:     # currently only 4 problem categories
                break
        fS.close()
        addToDProb(dNotCnv, dProb, 'NotConv', tK[0], tK[1])
        addToDProb(dSumHH, dProb, 'SumHetHom', tK[0], tK[1])
        addToDProb(dSumAP, dProb, 'SumAProp', tK[0], tK[1])

def getMode(lVs, lFs):
    assert len(lVs) == len(lFs)
    assert len(lVs) > 0
    cMd, cMxF = lVs[0], lFs[0]
    lMd = [cMd]
    for k in range(1, len(lVs)):
        if lFs[k] > cMxF:
            cMxF = lFs[k]
            cMd = lVs[k]
            lMd = [cMd]
        elif lFs[k] == cMxF:
            cMd = lVs[k]
            if cMd not in lMd:
                lMd.append(cMd)
    return lMd

def weighted_avg_and_std(arrVals, arrWts):
    """
    Return the weighted average and standard deviation.
    arrVals, arrWts -- Numpy ndarrays with the same shape.
    """
    assert arrVals.shape == arrWts.shape
    average = np.average(arrVals, weights = arrWts)
    variance = np.average((arrVals - average)**2, weights = arrWts)
    return (average, math.sqrt(variance))

def getStatVals(lVals, lFreqs):
    lMode = getMode(lVals, lFreqs)
    (wtMean, wtSD) = weighted_avg_and_std(np.array(lVals), np.array(lFreqs))
    return (lMode, wtMean, wtSD)

def getRoundPrec(valToRd, rdP):
    while len(str(round(valToRd, rdP))) >= 8:    # 8: length of tab
        rdP -= 1
    return rdP

def writeSumInfo(dInp, dInf, tD, nFT, rndP):
    createDir(os.path.join(tD, dInp['nD_ResultsC']))
    nFD = dInp['nF_Dic'] + '_' + dInp['calcMd']
    nFS = dInp['nF_Sum'] + '_' + dInp['calcMd']
    nFT += dInp['nFE_Txt']
    fTDic = open(os.path.join(tD, dInp['nD_ResultsC'], nFD + '_' + nFT), 'w')
    for cK, cV in sorted(dInf.items()):
        fTDic.write(str(cK) + ': ' + str(cV) + '\n')
    fTDic.close()
    fTSum = open(os.path.join(tD, dInp['nD_ResultsC'], nFS + '_' + nFT), 'w')
    dStatVs = {}
    for cMd in dInp['lNModel']:
        dStatVs[cMd] = {}
        fTSum.write('='*43 + ' ' + str(cMd) + ' ' + '='*43 + '\n')
        fTSum.write('-'*4 + ' Key: (max. intr. merit, min. intr. merit, box ' +
                    'seq. length mult., num. alleles init.) ' + '-'*4 + '\n')
        fTSum.write('Key\t\t\t' + 'Weighted mean\t' + 'Weighted SD\t' +
                    'Maximum\t\t' + 'Minimum\t\t' + 'Mode(s)\n')
        for setK in sorted(dInf['lVals'][cMd]):
            lVs, lFs = dInf['lVals'][cMd][setK], dInf['lFreqs'][cMd][setK]
            (lMode, wtMean, wtSD) = getStatVals(lVs, lFs)
            lMode = [round(cMode, rndP) for cMode in lMode]
            rPMn, rPSD = getRoundPrec(wtMean, rndP), getRoundPrec(wtSD, rndP)
            outSVs = [round(wtMean, rPMn), round(wtSD, rPSD),
                      round(dInf['maxV'][cMd][setK], rPMn),
                      round(dInf['minV'][cMd][setK], rPMn), lMode]
            dStatVs[cMd][setK] = outSVs
            fTSum.write(str(setK) + ':\t')
            for cSV in outSVs[:-1]:
                fTSum.write(str(cSV) + '\t\t')
            fTSum.write(str(outSVs[-1]) + '\n')
    fTSum.write('='*96 + '\n')
    fTSum.close()
    return dStatVs

def SUB_writeOtherInfo(dInp, dInf, fI, dSVs, keyInf, rdP, tpStr = 'Means'):
    iLInf = 0
    if tpStr == 'SDs':
        iLInf = 1
    fI.write('\n' + '*'*78 + ' ' + tpStr + ' ' + '*'*78 + '\n')
    dSVs[tpStr] = {}
    for cMd in dInp['lNModel']:
        dSVs[tpStr][cMd] = {}
        fI.write('='*78 + ' ' + str(cMd) + ' ' + '='*78 + '\n')
        fI.write('-'*39 + ' Key: (max. intr. merit, min. intr. merit,' +
                 'box seq. length mult., num. alleles init.) ' + '-'*39 + '\n')
        fI.write('Key\t\t\t' + 'Prop.Homozygotes\t' + 'Prop.Heterozygotes\t' +
                 'Av.Fit.Homozygotes\t' + 'Av.Fit.Heterozygotes\t' +
                 'HeterozygoteAdv.Abs.\t' + 'HeterozygoteAdv.Rel.\n')
        for setK in sorted(dInf[keyInf][cMd]):
            dVInf = dInf[keyInf][cMd][setK]
            lSVs = [dVInf[k][iLInf] for k in dVInf]
            for cI, cV in enumerate(lSVs):
                if cV != 'NA':
                    lSVs[cI] = round(cV, getRoundPrec(cV, rdP))
            dSVs[tpStr][cMd][setK] = lSVs
            fI.write(str(setK) + ':\t')
            for cSV in lSVs[:-1]:
                fI.write(str(cSV) + '\t\t\t')
            fI.write(str(lSVs[-1]) + '\n')

def writeOtherInfo(dInp, dInf, tD, nFT, rndP):
    createDir(os.path.join(tD, dInp['nD_ResultsC']))
    nFS, nFT = dInp['nF_Sum'] + '_' + dInp['calcMd'], nFT + dInp['nFE_Txt']
    fOInf = open(os.path.join(tD, dInp['nD_ResultsC'], nFS + '_' + nFT), 'w')
    dStatVs = {}
    SUB_writeOtherInfo(dInp, dInf, fOInf, dStatVs, 'HetAdvData', rndP, 'Means')
    SUB_writeOtherInfo(dInp, dInf, fOInf, dStatVs, 'HetAdvData', rndP, 'SDs')
    fOInf.close()
    return dStatVs

def writeProblems(dInp, dPrs, tD, nFT):
    createDir(os.path.join(tD, dInp['nD_ResultsC']))
    nFT += ('_' + dInp['calcMd'] + dInp['nFE_Txt'])
    fPrs = open(os.path.join(tD, dInp['nD_ResultsC'], nFT), 'w')
    for kL1 in ['NotConv', 'SumHetHom', 'SumAProp']:
        fPrs.write('*'*40 + ' ' + str(kL1) + ' ' + '*'*40 + '\n')
        fPrs.write('-'*4 + ' Key: (max. intr. merit, min. intr. merit, box ' +
                   'seq. length mult., num. alleles init.) ' + '-'*4 + '\n')
        for kL2 in sorted(dPrs[kL1]):
            fPrs.write('='*40 + ' ' + str(kL2) + ' ' + '='*40 + '\n')
            for kL3 in sorted(dPrs[kL1][kL2]):
                if kL1 == 'NotConv':
                    sVals = str(sorted(dPrs[kL1][kL2][kL3]))
                    fPrs.write(str(kL3) + ':\t' + sVals + '\n')
                else:
                    fPrs.write(str(kL3) + ':\n')
                    for cK, cV in sorted(dPrs[kL1][kL2][kL3].items()):
                        fPrs.write('\t' + str(cK) + ':\t' + str(cV) + '\n')
    fPrs.write('*'*90)
    fPrs.close()

def getCDiffCIn(v1, v2, lSmr, rP):
    lSmr[0] += 1
    cDf, cIn = 'NA', ''
    if v1 != 'NA' and v2 != 'NA':
        if v1 > 0.0:
            cDf = (v2 - v1)/v1*100.0
            cDf = round(cDf, getRoundPrec(cDf, rP))
        if v1 > v2:
            lSmr[1] += 1
            cIn = '+'
        elif v1 < v2:
            lSmr[2] += 1
            cIn = '-'
        else:
            lSmr[3] += 1
            cIn = '='
    return cDf, cIn

def removeNAs(dV):
    lK, lVN = [], []
    for cK, cV in dV.items():
        if cV != 'NA':
            lK.append(cK)
            lVN.append(cV)
    return (lK, lVN)

def plotModelComp(dInp, dDfs, sCt, tD, nFP):
    lPKeyF = ['rangeAFit', 'maxAFit', 'minAFit', 'multNumPatho', 'numAllIni']
    lPKeyL = ['allele intrinsic merit range', 'max. allele intrinsic merit',
              'min. allele intrinsic merit', 'length of "boxes" sequence',
              'initial number of alleles']
    lClrsDat = [(1.0, 0.0, 0.0), (0.9, 0.5, 0.0), (0.9, 0.9, 0.0),
                (0.0, 0.9, 0.0), (0.0, 0.2, 0.3), (0.0, 0.0, 0.9),
                (0.3, 0.0, 0.7), (0.6, 0.0, 0.4), (0.1, 0.1, 0.1)]
    (lKeys, yDt) = removeNAs(dDfs[sCt])
    nD_Plots = dInp['nD_Plots'] + '_' + dInp['calcMd']
    createDir(os.path.join(tD, nD_Plots))
    for k in range(len(lPKeyF)):
        if k == 0:
            xDt = [cK[k] - cK[k + 1] for cK in lKeys]
        else:
            xDt = [cK[k - 1] for cK in lKeys]
        nFPT = nFP + '_' + sCt + '_' + lPKeyF[k] + dInp['nFE_Plt']
        plt.figure(nFPT)
        plt.plot(xDt, yDt, color = lClrsDat[k], ls = 'None', marker = 'o',
                 label = '')
        plt.plot([min(xDt), max(xDt)], [0, 0], color = lClrsDat[-1], ls = '-')
        # a regression line might help
        fitFct = np.poly1d(np.polyfit(xDt, yDt, 1)) 
        plt.plot(xDt, fitFct(xDt), color = lClrsDat[-2], ls = '-')
        plt.title(nFP)
#         plt.legend(loc = 'best')
#         plt.ylim(ymin = 0)
        plt.xlabel(lPKeyL[k])
        plt.ylabel('model difference [' + sCt + '] (% dev. from ' +
                   dInp['lNModel'][0] + ' model)')
        plt.savefig(os.path.join(tD, nD_Plots, nFPT))
        plt.close()

def avDff(lDs, rdP):
    sMnDf = calcMeanNoNA(lDs)
    if sMnDf != 'NA':
        sMnDf = str(round(sMnDf, getRoundPrec(sMnDf, rdP)))
    return sMnDf

def writePlotModCmp(dInp, dS, tD, nFC, rndP, strCmp, iDt = -1):
    lNMs = dInp['lNModel']                 # for convenience (len(lNMs) >= 2)
    lIDt, lCats = [iDt], ['mean']
    if iDt < 0:
        lIDt, lCats = [0, 2, 3], ['mean', 'max', 'min']   # indices of lCats
    createDir(os.path.join(tD, dInp['nD_ResultsC']))
    fTCmp = open(os.path.join(tD, dInp['nD_ResultsC'], nFC), 'w')
    dDfs, plusL, dashL = {}, ('+'*90 + '\n'), ('-'*90 + '\n')
    fTCmp.write('='*22 + ' ' + strCmp + ' - model comparison ' + '='*22 + '\n')
    for iCat, sCat in enumerate(lCats):
        fTCmp.write(plusL)
        fTCmp.write('\t' + 'Rel. deviation of ' + sCat + ' of ' + lNMs[1] +
                    ' model to ' + sCat + ' of ' + lNMs[0] + ' model\n')
        fTCmp.write(plusL)
        fTCmp.write('-- Key: (max. intr. merit, min. intr. merit,' +
                    'box seq. length mult., num. alleles init.) -\n')
        fTCmp.write(dashL)
        fTCmp.write('Key\t\t\t' + 'Difference ' + sCat + ' value (rel./%)\t' +
                    sCat + ' value of ' + lNMs[0] + ' model higher?\n')
        dDfs[sCat] = {}
        lDfs, lTPME = [], [0, 0, 0, 0]
        for setK in sorted(dS[lNMs[0]]):
            if setK in dS[lNMs[1]]:
                iD = lIDt[iCat]
                vNRAsOD, vRAsOD = dS[lNMs[0]][setK][iD], dS[lNMs[1]][setK][iD]
                cDf, cIn = getCDiffCIn(vNRAsOD, vRAsOD, lTPME, rndP)
                dDfs[sCat][setK] = cDf
                lDfs.append(cDf)
                fTCmp.write(str(setK) + ':\t' + str(cDf) + '\t'*4 + cIn + '\n')
        fTCmp.write(dashL)
        fTCmp.write('Average difference of ' + sCat + ' (rel./%): ' +
                    avDff(lDfs, rndP) + '\n')
        fTCmp.write(dashL)
        fTCmp.write('Cases where ' + sCat + ' value of ' + lNMs[0] + ' model' +
                    ' is (greater, smaller, equal) to ' + sCat + ' value of ' +
                    lNMs[1] + ' model\n')
        fTCmp.write('(' + str(lTPME[1]) + ', ' + str(lTPME[2]) + ', ' +
                    str(lTPME[3]) + ')\n')
        fTCmp.write(dashL)
        fTCmp.write('Total cases: ' + str(lTPME[0]) + '\n')
        fTCmp.write(dashL)
        if dInp['dFlags']['plotMC']:
            plotModelComp(dInp, dDfs, sCat, tD, nFC[:-4])
    fTCmp.close()

def writeModelDistrComp(dInp, tD, nFC, rndP, strCmp):
    lNMs = dInp['lNModel']                 # for convenience (len(lNMs) >= 2)
    createDir(os.path.join(tD, dInp['nD_ResultsC']))
    fMDC = open(os.path.join(tD, dInp['nD_ResultsC'], nFC), 'w')
    fMDC.write('='*22 + ' ' + strCmp + ' - model distribution comparison ' +
               '='*22 + '\n')
    fMDC.close()

def extractAllInfo(dInp, tDir, rndPr):
    indDSRs = dInp['nD_ResultsC'] + '_' + dInp['calcMd'] + '_'
    # fill the dictionaries containing all summary information
    lInfCat = ['lVals', 'lFreqs', 'wMean', 'maxV', 'minV']
    dInfNA, dInfDF = {cK: {} for cK in lInfCat}, {cK: {} for cK in lInfCat}
    dOData = {'HetAdvData': {}}
    dProbs = {'NotConv': {}, 'SumHetHom': {}, 'SumAProp': {}}
    lDirs = os.listdir(tDir)
    nF_BSR = dInp['nF_BSR'] + dInp['nFE_Txt']
    nF_NAS = dInp['nF_NAS'] + dInp['nFE_Txt']
    nF_DFS = dInp['nF_DFS'] + dInp['nFE_Txt']
    nF_HAC = dInp['nF_HAC'] + dInp['nFE_Txt']
    nF_RPR = dInp['nF_RepRes'] + dInp['nFE_Txt']
    for cDir in lDirs:
        if len(cDir) >= len(indDSRs):
            if cDir[:len(indDSRs)] == indDSRs:
                extractDistributions(dInp, tDir, cDir, nF_BSR, nF_NAS, nF_DFS)
                extractValsFreqs(dInfNA, tDir, cDir, nF_NAS)
                extractValsFreqs(dInfDF, tDir, cDir, nF_DFS, False)
                extractOtherData(dOData, tDir, cDir, nF_HAC, dInp['dNumRep'])
                extractProblems(dProbs, tDir, cDir, nF_RPR)
    # extract useful information and write the dictionaries and summaries
    dSVNA = writeSumInfo(dInp, dInfNA, tDir, dInp['nF_NAT'], rndPr)
    dSVDF = writeSumInfo(dInp, dInfDF, tDir, dInp['nF_DFT'], rndPr)
    dSVOI = writeOtherInfo(dInp, dOData, tDir, dInp['nF_HAC'], rndPr)
    writeProblems(dInp, dProbs, tDir, dInp['nF_Prb'])
    nF_MCP = dInp['nF_MCP'] + '_' + dInp['calcMd']
    nF_MDC = dInp['nF_MDC'] + '_' + dInp['calcMd']
    if len(dInp['lNModel']) >= 2:
        nF1 = nF_MCP + '_' + dInp['nF_NAT'] + dInp['nFE_Txt']
        writePlotModCmp(dInp, dSVNA, tDir, nF1, rndPr, 'Number of alleles')
        nF2 = nF_MCP + '_' + dInp['nF_DFT'] + dInp['nFE_Txt']
        writePlotModCmp(dInp, dSVDF, tDir, nF2, rndPr,
                        'Intrinsic merit difference')
        for cI, cNCol in enumerate(dInp['lNF_HACol']):
            nFC = nF_MCP + '_' + cNCol + dInp['nFE_Txt']
            writePlotModCmp(dInp, dSVOI['Means'], tDir, nFC, rndPr, cNCol, cI)
        nF_MDC_1 = nF_MDC + '_' + dInp['nF_NAT'] + dInp['nFE_Txt']
        writeModelDistrComp(dInp, tDir, nF_MDC_1, rndPr, 'Number of alleles')
    print('Finished extracting info.')

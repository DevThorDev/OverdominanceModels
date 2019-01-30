################################################################################
# WriteAndPrint.py #
################################################################################

import os, pprint

from Constants import C_TSTEPPING, RPR_15, RPR_12, RPR_04
from SpecialFunctions import getPRes, updateResDict, getFittestAllele
from GenericFunctions import (createDir, calcRunMean, addEntryToDict,
                              addEntryToDVList)

def deleteAllFiles(dictI, cDir):
    for cMdl in dictI['lModels']:
        pRes = cDir + '/' + dictI['nD_ResultsC'] + getPRes(dictI, cMdl)
        createDir(pRes)
        for nOneF in os.listdir(pRes):
            os.remove(pRes + '/' + nOneF)

def printCurTask(iT, nT, cDT, nMd = ''):
    strPr = '-'*8 + ' Model: ' + str(nMd) + ' - task:'
    if len(nMd) == 0:
        print('*'*32, 'Starting task', iT + 1, 'of', nT, '*'*32)
        strPr = '-'*16 + ' Task:'
    for cK, cV in sorted(cDT.items()):
        strPr += (' ' + str(cK) + '_' + str(cV))
    strPr += (' ' + '-'*8)
    if len(nMd) == 0:
        strPr += ('-'*(16 - 8))
    print(strPr)

def printCurRepet(dictI, iM, iT, nM, nT, cRep):
    percRep = round((cRep - 1)*100.0/dictI['numRep'], 2)
    percTotal = round(iT*100.0/nT + (cRep - 1)*100.0/(nT*dictI['numRep']) +
                      iM*100.0/(nT*dictI['numRep']*nM), 2)
    print('+'*12 + ' Repetition', cRep, 'of', dictI['numRep'], '+'*12 + ' (',
          percRep, '% r. / ', percTotal, '% t.) ' + '+'*12)

def printTaskFinished(iT, nT):
    print('*'*20, 'Finished task', iT + 1, 'of', nT, '(' +
          str(round((iT + 1)*100.0/nT, 2)), '% of tasks done)', '*'*20)

def openOutFiles(dictI, dictGP, cRp, cDir):
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    # open step result file
    fRes = open(pRes + '/' +
                dictI['nF_StepRes'] + '_Rep' + str(cRp) + dictI['nFE_Txt'], 'w')
    fRes.write('TimeStep PopulationFitness PropHomozygotes PropHeterozygotes')
    for curAll in sorted(dictGP):
        fRes.write(' PropAllele' + str(curAll))
    fRes.write('\n')
    fRes.close()
    # open rate result file
    fRes = open(pRes + '/' +
                dictI['nF_RateRes'] + '_Rep' + str(cRp) + dictI['nFE_Txt'], 'w')
    fRes.write('TimeStep NumberAllelesCurrent NumberAllelesLostTimeStep' +
               ' NumberAllelesLostTotal\n')
    fRes.close()
    
def writeToStepResult(dictI, dictGP, dictR, cRp, tSt, cDir):
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    fRes = open(pRes + '/' +
                dictI['nF_StepRes'] + '_Rep' + str(cRp) + dictI['nFE_Txt'], 'a')
    fRes.write(str(tSt) + ' ' + str(dictR['popFitness']))
    fRes.write(' ' + str(dictR['propZyg']['propHomoZyg']) +
               ' ' + str(dictR['propZyg']['propHeteroZyg']))
    for curAll in range(dictI['numAllIni']):
        try:
            fRes.write(' ' + str(dictGP[curAll][2]))
        except KeyError:
            fRes.write(' 0.0')
    fRes.write('\n')
    fRes.close()

def writeToRateResult(dictI, dictR, cRp, tSt, cDir):
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    fRes = open(pRes + '/' +
                dictI['nF_RateRes'] + '_Rep' + str(cRp) + dictI['nFE_Txt'], 'a')
    fRes.write(str(tSt) + ' ' + str(len(dictR['lAllPers'])) +
               ' ' + str(dictR['nAllTSLost']) +
               ' ' + str(dictR['nAllTotLost']) + '\n')
    fRes.close()

def printGPResult(dictGP, dictR, tSt):
    print('--------------------------------------------------')
    print('Time step', tSt, '- Population Fitness', dictR['popFitness'],
          '- Genotype dictionary:')
    pprint.pprint(dictGP)

def writeGPResult(dictI, dictGP, cRp, tSt, cDir):
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    fGPRes = open(pRes + '/' +
                  dictI['nF_GPRes'] + '_Rep' + str(cRp) + '_timeStep' +
                  str(tSt) + dictI['nFE_Txt'], 'w')
    fGPRes.write('Allele AlleleCharacteristic AlleleFitness AlleleProportion' +
                 ' MarginalFitness\n')
    for dEntry in dictGP.items():
        fGPRes.write(str(dEntry[0]) + ' ')
        for oneEl in dEntry[-1][:-1]:
            fGPRes.write(str(oneEl) + ' ')
        fGPRes.write(str(dEntry[-1][-1]) + '\n')
    fGPRes.close()

def writeGPShort(dictI, dictGP, dictR, cRp, tSt, cDir):
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    fSGPRes = open(pRes + '/' +
                   dictI['nF_SGPRes'] + '_Rep' + str(cRp) + '_timeStep' +
                   str(tSt) + dictI['nFE_Txt'], 'w')
    fSGPRes.write('Count Allele' +
                  ' AlleleFitness' +
                  ' AlleleProportion MarginalFitness\n')
    lAllUsed = []
    updateResDict(dictI, dictGP, dictR)
    for persAllIdx in range(len(dictR['lAllPers'])):
        allMx = getFittestAllele(dictGP, dictR['lAllPers'], lAllUsed)
        fSGPRes.write(str(persAllIdx + 1) + ' ' + str(allMx) + ' ')
        for oneEl in dictGP[allMx][1:-1]:
            fSGPRes.write(str(oneEl) + ' ')
        fSGPRes.write(str(dictGP[allMx][-1]) + '\n')
    fSGPRes.close()

def getHetAdvStr(dHA, keyHom, keyHet):
    sHAAbs, sHARel = 'NA', 'NA'
    if dHA[keyHet] != 'NA':
        sHAAbs = str(dHA[keyHet] - dHA[keyHom])
        if dHA[keyHom] > 0:
            sHARel = str((dHA[keyHet] - dHA[keyHom])/dHA[keyHom]*100.0)
    return sHAAbs, sHARel

def writeFinalData(dictI, dictGP, dictGTF, dictRR, cRp, cDir):
    dFin_S = {'numZygotes': {'numHom': 0, 'numHet': 0},
              'hetAdv': {'avFHom': 0.0, 'avFHet': 0.0, 'curHetAdv': 0.0},
              'hetOverlp': {'sumHetOver': 0.0, 'avHetOver': 0.0}}
    dFin_E = {'numZygotes': {'numHom': 0, 'numHet': 0},
              'hetAdv': {'avFHom': 0.0, 'avFHet': 0.0, 'curHetAdv': 0.0},
              'hetOverlp': {'sumHetOver': 0.0, 'avHetOver': 0.0}}
    dFin_C = {'propOfHom': 0.0, 'propOfHet': 0.0, 'avFitHom': 0.0,
              'avFitHet': 0.0}
    for oneKey, oneVal in dictGTF.items():
        if oneKey[0] == oneKey[1]:  # homozygote
            dFin_S['numZygotes']['numHom'] += 1
            avFitHom = calcRunMean(oneVal[0], dFin_S['hetAdv']['avFHom'],
                                   dFin_S['numZygotes']['numHom'])
            dFin_S['hetAdv']['avFHom'] = avFitHom
            if oneKey[0] in dictGP and oneKey[1] in dictGP:
                dFin_E['numZygotes']['numHom'] += 1
                avFitHom = calcRunMean(oneVal[0], dFin_E['hetAdv']['avFHom'],
                                       dFin_E['numZygotes']['numHom'])
                dFin_E['hetAdv']['avFHom'] = avFitHom
        else:                       # heterozygote
            dFin_S['numZygotes']['numHet'] += 1
            avFitHet = calcRunMean(oneVal[0], dFin_S['hetAdv']['avFHet'],
                                   dFin_S['numZygotes']['numHet'])
            dFin_S['hetAdv']['avFHet']= avFitHet
            dFin_S['hetOverlp']['sumHetOver'] += oneVal[1]
            if oneKey[0] in dictGP and oneKey[1] in dictGP:
                dFin_E['numZygotes']['numHet'] += 1
                avFitHet = calcRunMean(oneVal[0], dFin_E['hetAdv']['avFHet'],
                                       dFin_E['numZygotes']['numHet'])
                dFin_E['hetAdv']['avFHet']= avFitHet
                dFin_E['hetOverlp']['sumHetOver'] += oneVal[1]
    startHetAdv = dFin_S['hetAdv']['avFHet'] - dFin_S['hetAdv']['avFHom']
    dFin_S['hetAdv']['curHetAdv'] = startHetAdv
    endHetAdv = dFin_E['hetAdv']['avFHet'] - dFin_E['hetAdv']['avFHom']
    dFin_E['hetAdv']['curHetAdv'] = endHetAdv
    if (dFin_S['numZygotes']['numHet']) > 0:
        dFin_S['hetOverlp']['avHetOver'] = (dFin_S['hetOverlp']['sumHetOver']/
                                             dFin_S['numZygotes']['numHet'])
    else:
        dFin_S['hetOverlp']['avHetOver'] = 0.0
    if (dFin_E['numZygotes']['numHet']) > 0:
        dFin_E['hetOverlp']['avHetOver'] = (dFin_E['hetOverlp']['sumHetOver']/
                                             dFin_E['numZygotes']['numHet'])
    else:
        dFin_E['hetOverlp']['avHetOver'] = 0.0
    for cAllID in dictGP:
        for oAllID in dictGP:
            propGT = dictGP[cAllID][2]*dictGP[oAllID][2]
            keyGTF = (max(cAllID, oAllID), min(cAllID, oAllID))
            if keyGTF in dictGTF:           # otherwise disregard this allele
                if oAllID == cAllID:        # homozygote
                    dFin_C['propOfHom'] += propGT
                    dFin_C['avFitHom'] += propGT*dictGTF[keyGTF][0]
                else:                       # heterozygote
                    dFin_C['propOfHet'] += propGT
                    dFin_C['avFitHet'] += propGT*dictGTF[keyGTF][0]
    if dFin_C['propOfHom'] > 0.0:
        dFin_C['avFitHom'] /= dFin_C['propOfHom']
    else:
        dFin_C['avFitHom'] = 'NA'
    if dFin_C['propOfHet'] > 0.0:
        dFin_C['avFitHet'] /= dFin_C['propOfHet']
    else:
        dFin_C['avFitHet'] = 'NA'
    cNMod = dictI['dNModel'][dictI['domType']]
    dictRR[cNMod]['dHetAdvCorrect'][cRp] = dFin_C
    sumPropHHRd = round(dFin_C['propOfHom'] + dFin_C['propOfHet'], RPR_12)
    if sumPropHHRd != 1.0:
        print('WARNING: Proportions of homozygotes and heterozygotes sum to',
              sumPropHHRd)
        dictRR[cNMod]['Problems']['SumHetHom'][cRp] = sumPropHHRd
    pRes = cDir + '/' + dictI['nD_Results']
    createDir(pRes)
    if dictI['dFlags']['saveAll']:
        fHARes = open(pRes + '/' + dictI['nF_HetAdv'] +
                      '_Rep' + str(cRp) + dictI['nFE_Txt'], 'w')
        fHARes.write('StartNumHomozygotes StartNumHeterozygotes' +
                     ' StartAvFitnessHomozygotes StartAvFitnessHeterozygotes' +
                     ' StartHeterozygoteAdvantage StartAvHeterozygoteOverlap' +
                     ' EndNumHomozygotes EndNumHeterozygotes' +
                     ' EndAvFitnessHomozygotes EndAvFitnessHeterozygotes' +
                     ' EndHeterozygoteAdvantage EndAvHeterozygoteOverlap\n')
        fHARes.write(str(dFin_S['numZygotes']['numHom']) + ' ' +
                     str(dFin_S['numZygotes']['numHet']) + ' ' +
                     str(dFin_S['hetAdv']['avFHom']) + ' ' +
                     str(dFin_S['hetAdv']['avFHet']) + ' ' +
                     str(startHetAdv) + ' ' +
                     str(dFin_S['hetOverlp']['avHetOver']) + ' ' +
                     str(dFin_E['numZygotes']['numHom']) + ' ' +
                     str(dFin_E['numZygotes']['numHet']) + ' ' +
                     str(dFin_E['hetAdv']['avFHom']) + ' ' +
                     str(dFin_E['hetAdv']['avFHet']) + ' ' +
                     str(endHetAdv) + ' ' +
                     str(dFin_E['hetOverlp']['avHetOver']) + '\n')
        fHARes.close()
        fHARes = open(pRes + '/' + dictI['nF_HAC'] +
                      '_Rep' + str(cRp) + dictI['nFE_Txt'], 'w')
        for cNCol in dictI['lNF_HACol'][:-1]:
            fHARes.write(cNCol + ' ')
        fHARes.write(dictI['lNF_HACol'][-1] + '\n')
        sHetAdvAbs, sHetAdvRel = getHetAdvStr(dFin_C, 'avFitHom', 'avFitHet')
        if dFin_C['avFitHet'] != 'NA':
            sHetAdvAbs = str(dFin_C['avFitHet'] - dFin_C['avFitHom'])
            if dFin_C['avFitHom'] > 0:
                sHetAdvRel = str((dFin_C['avFitHet'] - dFin_C['avFitHom'])/
                                   dFin_C['avFitHom']*100.0)
        fHARes.write(str(dFin_C['propOfHom']) + ' ' +
                     str(dFin_C['propOfHet']) + ' ' +
                     str(dFin_C['avFitHom']) + ' ' +
                     str(dFin_C['avFitHet']) + ' ' +
                     sHetAdvAbs + ' ' + sHetAdvRel + '\n')
        fHARes.close()

def removeCurRepFiles(dictI, cRp, pRes):
    for nOneF in os.listdir(pRes):
        lFBits = nOneF.split('.')[0].split('_')
        for oneBit in lFBits:
            if oneBit[:3] == 'Rep':
                if oneBit[3:] == str(cRp):
                    if os.path.isfile(pRes + '/' + nOneF):
                        os.remove(pRes + '/' + nOneF)

def checkNumAllPers(dI, dR, dRR, cRp, cMx = 0):
    cNMod = dI['dNModel'][dI['domType']]
    print('Number of alleles persisting:', dR['nAllPers'], '(previous max.:',
          str(cMx) + ')')
    if dR['nAllPers'] >= dI['numAllThr'] + dI['dClOutc']['huge']:
        dRR[cNMod]['clOutcInc']['huge'] += 1
    elif dR['nAllPers'] >= dI['numAllThr'] + dI['dClOutc']['large']:
        dRR[cNMod]['clOutcInc']['large'] += 1
    elif dR['nAllPers'] <= dI['numAllThr'] + dI['dClOutc']['wee']:
        dRR[cNMod]['clOutcInc']['wee'] += 1
    elif dR['nAllPers'] <= dI['numAllThr'] + dI['dClOutc']['small']:
        dRR[cNMod]['clOutcInc']['small'] += 1
    else:
        dRR[cNMod]['clOutcInc']['normal'] += 1
    dRR[cNMod]['dAllPers'][cRp] = dR['nAllPers']

def writeResults(dictI, dictGP, dictGTF, dictR, dictRR, curRp, tSt, cDir, iC):
    fGoOn, cMNA, cNMod = True, 0, dictI['dNModel'][dictI['domType']]
    if tSt > 0 and len(dictRR[cNMod]['dAllPers']) > 0:
        cMNA = max(list(dictRR[cNMod]['dAllPers'].values()))
    if not iC:
        if ((tSt == 0 or (tSt > 0 and dictR['nAllPers'] >= cMNA) or
             not dictI['dFlags']['breakEarly']) and dictR['nAllPers'] > 1):
            if (tSt%dictI['modUpdt'] == 0 or tSt%dictI['modStep'] == 0 or
                tSt%dictI['modRate'] == 0):
                updateResDict(dictI, dictGP, dictR)
            if tSt%dictI['modRate'] == 1:
                dictR['nAllTSLost'] = 0
            if tSt%dictI['modDisp'] == 0 and tSt > 0:
                print('Repetition', curRp, 'time step', tSt, '- number of',
                      'alleles remaining:', dictR['nAllPers'], '(not conv. /',
                      'current max.:', str(cMNA) + ')')
            if tSt%dictI['modStep'] == 0 and dictI['dFlags']['saveAll']:
                writeToStepResult(dictI, dictGP, dictR, curRp, tSt, cDir)
            if tSt%dictI['modRate'] == 0 and dictI['dFlags']['saveAll']:
                writeToRateResult(dictI, dictR, curRp, tSt, cDir)
            if tSt%dictI['modFile'] == 0 and dictI['dFlags']['saveAll']:
                writeGPResult(dictI, dictGP, curRp, tSt, cDir)
                writeGPShort(dictI, dictGP, dictR, curRp, tSt, cDir)
            if tSt == dictI['numIter']:
                writeFinalData(dictI, dictGP, dictGTF, dictRR, curRp, cDir)
                checkNumAllPers(dictI, dictR, dictRR, curRp, cMNA)
        else:
            fGoOn = False
            removeCurRepFiles(dictI, curRp, cDir + '/' + dictI['nD_Results'])
            sPr = ('Stopping repetition ' + str(curRp) + ' at time step ' +
                   str(tSt) + ' as current number of alleles = ' +
                   str(dictR['nAllPers']))
            if dictR['nAllPers'] < cMNA:
                sPr += (' < ' + str(cMNA) + ' (current max.)')
            print(sPr)
        if tSt == dictI['numIter']:
            dictRR[cNMod]['Problems']['NotConv'].append(curRp)
    else:
        fGoOn = False
        if dictI['calcMd'] == C_TSTEPPING:
            print('Converged at time step', tSt, 'of', dictI['numIter'], '-',
                  'number of alleles remaining:', dictR['nAllPers'])
        dictR['nAllTSLost'] = 0
        updateResDict(dictI, dictGP, dictR)
        if dictI['dFlags']['saveAll']:
            writeToStepResult(dictI, dictGP, dictR, curRp, tSt, cDir)
            writeToRateResult(dictI, dictR, curRp, tSt, cDir)
            writeGPResult(dictI, dictGP, curRp, tSt, cDir)
            writeGPShort(dictI, dictGP, dictR, curRp, tSt, cDir)
        writeFinalData(dictI, dictGP, dictGTF, dictRR, curRp, cDir)
        checkNumAllPers(dictI, dictR, dictRR, curRp, cMNA)
        if (cMNA > 0 and dictR['nAllPers'] >= cMNA) or curRp == 1:
            addEntryToDVList(dictRR[cNMod]['dMaxNA'], dictR['nAllPers'], curRp)
    return fGoOn

def writeDistsRes(dI, dDistNA, dDistDF, cDir):
    for cMdl in dI['lModels']:
        pRes = cDir + '/' + dI['nD_ResultsC'] + getPRes(dI, cMdl)
        createDir(pRes)
        # open number of alleles distribution file
        dDistNAM = dDistNA[dI['dNModel'][cMdl]]
        fRes = open(pRes + '/' + 'DistrNumAllelesResult' + dI['nFE_Txt'], 'w')
        fRes.write('NumberAlleles TimesRecorded\n')
        for cKey in sorted(dDistNAM, reverse = True):
            fRes.write(str(cKey) + ' ' + str(dDistNAM[cKey]) + '\n')
        fRes.close()
        # open fitness difference distribution file
        dDistDFM = dDistDF[dI['dNModel'][cMdl]]
        fRes = open(pRes + '/' + 'DistrDiffFitResult' + dI['nFE_Txt'], 'w')
        fRes.write('FitnessDifferenceAbs TimesRecorded\n')
        for cKey in sorted(dDistDFM, reverse = True):
            fRes.write(str(round(cKey, RPR_12)) + ' ' + str(dDistDFM[cKey]) +
                       '\n')
        fRes.close()

def writeBasicRes(dI, dBasR, cDir):
    for cMdl in dI['lModels']:
        pRes = cDir + '/' + dI['nD_ResultsC'] + getPRes(dI, cMdl)
        createDir(pRes)
        # open basic result file
        dBasRM = dBasR[dI['dNModel'][cMdl]]
        fRes = open(pRes + '/' + 'BasicResult' + dI['nFE_Txt'], 'w')
        fRes.write('Repetition NumberAlleles FitnessDifferenceAbs\n')
        for cKey in sorted(dBasRM):
            fRes.write(str(cKey) + ' ' + str(dBasRM[cKey][0]) + ' ' +
                       str(dBasRM[cKey][1]) + '\n')
        fRes.close()

def writeHetAdvCorrectRes(dI, dRR, cDir):
    for cMdl in dI['lModels']:
        cNMod = dI['dNModel'][cMdl]
        dHACM = dRR[cNMod]['dHetAdvCorrect']
        # if breaking early, clean the heterozygote advantage dictionary
        if dI['dFlags']['breakEarly'] and len(dRR[cNMod]['dMaxNA']) > 0:
            for cRp in dRR[cNMod]['dAllPers']:
                if dRR[cNMod]['dAllPers'][cRp] < max(dRR[cNMod]['dMaxNA']):
                    if cRp in dRR[cNMod]['dHetAdvCorrect']:
                        del dRR[cNMod]['dHetAdvCorrect'][cRp]
        # create directory and open heterozygote advantage correct file
        pRes = cDir + '/' + dI['nD_ResultsC'] + getPRes(dI, cMdl)
        createDir(pRes)
        fRes = open(pRes + '/' + dI['nF_HAC'] + dI['nFE_Txt'], 'w')
        fRes.write('Repetition ProportionHomozygotes ProportionHeterozygotes' +
                   ' AvFitHomozygotes AvFitHeterozygotes' +
                   ' HeterozygoteAdvAbs HeterozygoteAdvRel\n')
        for cRep in sorted(dHACM):
            cDHA = dHACM[cRep]
            sHetAdvAbs, sHetAdvRel = getHetAdvStr(cDHA, 'avFitHom', 'avFitHet')
            fRes.write(str(cRep) + ' ')
            fRes.write(str(cDHA['propOfHom']) + ' ' +
                       str(cDHA['propOfHet']) + ' ' +
                       str(cDHA['avFitHom']) + ' ' +
                       str(cDHA['avFitHet']) + ' ' +
                       sHetAdvAbs + ' ' + sHetAdvRel + '\n')
        fRes.close()

def getRepMaxNumAllPers(dNAP):
    maxNAP = max(dNAP)
    lRepMaxNAP = dNAP[maxNAP]
    return maxNAP, lRepMaxNAP

def delAllNonHuge(dI, dRR, cMdl, pRes):
    dAllPersM, cMNA = dRR[dI['dNModel'][cMdl]]['dAllPers'], 0
    if len(dAllPersM) > 0:
        cMNA = max(list(dAllPersM.values()))
    for cRp in dAllPersM:
        if dAllPersM[cRp] < cMNA:
            removeCurRepFiles(dI, cRp, pRes)

def writeRepRes(dictI, dictRR, cDir):
    for cMdl in dictI['lModels']:
        dictRRM = dictRR[dictI['dNModel'][cMdl]]
        pRes = cDir + '/' + dictI['nD_ResultsC'] + getPRes(dictI, cMdl)
        createDir(pRes)
        fCur = open(pRes + '/' + dictI['nF_RepRes'] + dictI['nFE_Txt'], 'w')
        sCL80 = '-'*80 + '\n'
        fCur.write(sCL80)
        fCur.write('Repetitions that did not converge:\n')
        if len(dictRRM['Problems']['NotConv']) > 0:
            fCur.write(str(dictRRM['Problems']['NotConv']) + '\n')
        fCur.write(sCL80)
        fCur.write('Repetitions where frequencies of homo- and' +
                   ' heterozygotes did not add up to 1:\n')
        for cK, cV in sorted(dictRRM['Problems']['SumHetHom'].items()):
            fCur.write(str(cK) + ':\t' + str(cV) + '\n')
        fCur.write(sCL80)
        fCur.write('Repetitions where allele frequencies did not add up to' +
                   ' 1:\n')
        for cK, cV in sorted(dictRRM['Problems']['SumAProp'].items()):
            fCur.write(str(cK) + ':\t' + str(cV) + '\n')
        fCur.write(sCL80)
        lClOutc = sorted(dictI['dClOutc'])
        for oneKey in lClOutc[:-1]:
            fCur.write(oneKey + ' ')
        fCur.write(lClOutc[-1] + '\n')
        fCur.write(sCL80)
        for oneKey in lClOutc[:-1]:
            fCur.write(str(dictRRM['clOutcInc'][oneKey]) + ' ')
        fCur.write(str(dictRRM['clOutcInc'][lClOutc[-1]]) + '\n')
        fCur.write(sCL80)
        if len(dictRRM['dMaxNA']) > 0:
            fCur.write('List of repetition numbers for maximum number of' +
                       ' alleles persisting:\n')
            keyDMNAP, lValDMNAP = getRepMaxNumAllPers(dictRRM['dMaxNA'])
            fCur.write(str(keyDMNAP) + ':\t' + str(lValDMNAP) + '\n')
            fCur.write(sCL80)
        fCur.close()
        if dictI['dFlags']['delNonHuge']:
            delAllNonHuge(dictI, dictRR, cMdl, pRes)

def writeSglTaskVals(cD, cK, cL):
    strL = ''
    if cK == 'maxAFit':
        if 'minAFit' in cD:
            for cI in range(len(cL)):
                strL += (str(max(cD[cK][cI], cD['minAFit'][cI])) + ', ')
    elif cK == 'minAFit':
        if 'maxAFit' in cD:
            for cI in range(len(cL)):
                strL += (str(min(cD[cK][cI], cD['maxAFit'][cI])) + ', ')
    else:
        for cEl in cL:
            strL += (str(cEl) + ', ')
    return strL

def writeLTasks(dictI, lDT, cDir):
    dOut, nT = {}, 0
    for cDT in lDT:
        for cK, cV in sorted(cDT.items()):
            addEntryToDict(dOut, cK, cV)
    pRes = cDir + '/' + dictI['nD_ResultsC']
    createDir(pRes)
    fTL = open(pRes + '/' + dictI['nF_TKL'] + dictI['nFE_Txt'], 'w')
    for cKOut, cVOut in sorted(dOut.items()):
        nT = max(nT, len(cVOut))
        strOutTL = str(cKOut) + ' = c('
        strOutTL += writeSglTaskVals(dOut, cKOut, cVOut)
        strOutTL = strOutTL[:-2] + ')\n'
        fTL.write(strOutTL)
    fTL.close()
    fStSz = open(pRes + '/' + dictI['nF_STS'] + dictI['nFE_Txt'], 'w')
    strOutStSz = 'stepSizes = c('
    for k in range(nT):
        cStSz = round(abs(dOut['maxAFit'][k] - dOut['minAFit'][k])/
                      dOut['numAllIni'][k], RPR_15)
        strOutStSz += (str(cStSz) + ', ')
    strOutStSz = strOutStSz[:-2] + ')\n'
    fStSz.write(strOutStSz)
    fStSz.close()

def printElapsedTimeSim(stT, cT, sPre = 'Time'):
    # calculate and display elapsed time 
    elT = round(cT - stT, RPR_04)
    print(sPre, 'elapsed:', elT, 'seconds, this is', round(elT/60, RPR_04),
          'minutes or', round(elT/3600, RPR_04), 'hours or',
          round(elT/(3600*24), RPR_04), 'days.')

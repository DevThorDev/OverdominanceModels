################################################################################
# GenericFunctions.py #
################################################################################

import os, random, math

import numpy as np

def createDir(pF):
    if not os.path.isdir(pF):
        os.mkdir(pF)

def drawFromPList(lProb):
    if len(lProb) == 0:
        print('Error: Probability vector has size 0! Returning 0...')
        return 0
    drIdx = 0
    curSum = 0.0
    drNum = random.random()
    while curSum <= drNum and drIdx < len(lProb):
        curSum += lProb[drIdx]
        drIdx += 1
    return drIdx - 1

def addToCountDict(theDict, cEl, cInc = 1):
    if cEl in theDict:
        theDict[cEl] += cInc
    else:
        theDict[cEl] = cInc

def convSubSetToBinChain(lA, nE, simplfyIt = False):
    dA = {}
    # convert to chains of 1s and 0s
    for lIdx in range(len(lA)):
        lRespPat = [0]*nE
        n1s = 0
        for k in lA[lIdx]:
            lRespPat[k] = 1
            n1s += 1
        dA[lIdx] = [lRespPat, n1s/nE]
        if simplfyIt:
            dA[lIdx][0] = n1s/nE
    return dA

def randPartitionRange(numParts, rLen = 1.0):
    lShares = [0.0]*numParts
    for i in range(numParts):
        lShares[i] = random.random()
    sumShares = math.fsum(lShares)
    for i in range(numParts):
        lShares[i] *= rLen/sumShares
    return lShares

def calcRelFitness(l1, l2):
    lenGT = 0
    ovGT = 0
    for i in range(len(l1)):
        if l1[i] == 1 or l2[i] == 1:
            lenGT += 1
            if l1[i] == 1 and l2[i] == 1:
                ovGT += 1
    return [lenGT/len(l1), ovGT/len(l1)]

def calcMean(lNums):
    if len(lNums) > 0:
        return sum(lNums)/float(len(lNums))
    else:
        print('Error: list is empty - cannot calculate arithmetic mean.')
        return 0.0;

def calcRunMean(curVal, meanBef, numValsCur):
    if numValsCur > 0:
        return meanBef + (curVal - meanBef)/numValsCur
    else:
        return 0.0

def getKeyOfMaxDictVal(theDict, lKFull, lKUsed, dIdx, outStr):
    keyMaxVal = -1
    maxVal = 0.0
    for oneKey in lKFull:
        if oneKey not in lKUsed:
            if theDict[oneKey][dIdx] >= maxVal:
                keyMaxVal = oneKey
                maxVal = theDict[oneKey][dIdx]
    if keyMaxVal >= 0:
        lKUsed.append(keyMaxVal)
        return keyMaxVal
    else:
        print(outStr)
        return 0

def addEntryToDict(theDict, theKey, valToAdd):
    if theKey in theDict:
        theDict[theKey].append(valToAdd)
    else:
        theDict[theKey] = [valToAdd]

def addEntryToDVList(theDict, theKey, valToAdd):
    if theKey in theDict:
        valList = theDict[theKey]
        valList.append(valToAdd)
        theDict[theKey] = valList
    else:
        theDict[theKey] = [valToAdd]

def getNElD(inD):
    nElTot, nCt = 1, len(inD)
    lNElCt = [0]*nCt
    for (cI, (_, cV)) in enumerate(sorted(inD.items())):
        nElTot *= len(cV)
        lNElCt[cI] = len(cV)
    return nElTot, lNElCt, nCt

def addNewD(inD, nLD, lIs, nCt, iElT):
    nD = {}
    for (cI, (cK, cV)) in enumerate(sorted(inD.items())):
        nD[cK] = cV[lIs[cI]]
    nLD[iElT] = nD

def modLIdx(lIs, lNElCt, nCt):
    iCt = 0
    while lIs[iCt] >= lNElCt[iCt] - 1:
        if iCt < nCt - 1:
            iCt += 1
        else:
            iCt = 0
    lIs[iCt] += 1
    for iP in range(iCt):
        lIs[iP] = 0

def multDictVal(theDict, theKey, lIdx = -1, theMult = 1):
    if lIdx >= 0:
        theDict[theKey][lIdx] *= theMult
    else:
        theDict[theKey] *= theMult

def safeRemoveEl(theDict, theKey, lIdx = -1, elSzRem = 0.0):
    del theDict[theKey]
    if elSzRem != 0.0:
        cMult = 1.0
        if elSzRem < 0:
            cMult = 1.0/(1.0 + elSzRem)
        else:
            cMult = 1.0/(1.0 - elSzRem)
        for oneKey in theDict:
            multDictVal(theDict, oneKey, lIdx, cMult)

def removeElVect(vIn, lIRt = []):
    iVOut, vOut = 0, np.zeros(len(lIRt))
    for iVIn in range(vIn.shape[0]):
        if iVIn in lIRt:
            vOut[iVOut] = vIn[iVIn]
            iVOut += 1
    return vOut

def checkIfSqMat(cMat):
    isSqMat = False
    if cMat.ndim == 2:
        if cMat.shape[0] == cMat.shape[1]:
            isSqMat = True
    return isSqMat

def shrinkMat(matIn, lIRt = [], cAx = 0):
    assert cAx in [0, 1]
    iRCOut = 0
    matOut = np.zeros((len(lIRt), matIn.shape[0]))
    if cAx == 1:        # shrink by removing columns
        matOut = np.zeros((len(lIRt), len(lIRt)))
    for iRCIn in range(matIn.shape[cAx]):
        if iRCIn in lIRt:
            if cAx == 0:
                matOut[iRCOut, :] = matIn[iRCIn, :]
            else:
                matOut[:, iRCOut] = matIn[:, iRCIn]
            iRCOut += 1
    return matOut

def removeRowColMat(matIn, lIRt = []):
    assert matIn.ndim == 2
    assert matIn.shape[0] == matIn.shape[1]
    return shrinkMat(shrinkMat(matIn, lIRt), lIRt, 1)

def checkIfSuccess(lFlag, nElToCheck = 1, doAND = True):
    numSucc, nElToCheck = 0, min(nElToCheck, len(lFlag))
    for cFlag in lFlag[:nElToCheck]:
        if cFlag:
            numSucc += 1
    if doAND:
        return numSucc == nElToCheck
    else:
        return numSucc > 0

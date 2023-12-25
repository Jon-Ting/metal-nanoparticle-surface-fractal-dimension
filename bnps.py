from os import getcwd, listdir
from os.path import isdir
import pickle
import sys

from natsort import natsorted
import numpy as np
import pandas as pd

from sphractal.boxCnt import runBoxCnt


# Important variables
testCasePath = 'testCases/BNPs'
voxelSurf, exactSurf = False, True
numPoints, numBoxLen, minSample = 2500, 10, 5
radType, radMult, alphaMult, bulkCN = 'atomic', 1.2, 2.0, 12
minLenMult, maxLenMult = 0.15, 2.00

# General variables
trimLen, calcBL, confLvl = True, False, 95
vis, figType, rmInSurf, showPlot, writeBox, saveFig, verbose = True, 'paper', True, False, True, True, True
gridNum, findSurfAlg, genPCD = 1024, 'alphaShape', False
bufferDist, numCPUs = 5.0, int(sys.argv[1])

PROJECT_DIR = getcwd()
FASTBC = f"{PROJECT_DIR}/fastbc/3DbinImBCcpu.exe"
OUTPUT_DIR = f"{PROJECT_DIR}/outputs"


if __name__ == '__main__':

    if verbose: 
        print('Collecting test case xyz paths...')
    testCases = []
    for bnpName in listdir(testCasePath):
        for fName in listdir(f"{testCasePath}/{bnpName}"):
            if not isdir(f"{testCasePath}/{bnpName}/{fName}") and fName.endswith('xyz'):
                testCases.append(f"{testCasePath}/{bnpName}/{fName}")

    if verbose: 
        print('Running once for JIT compilation...')
    _ = runBoxCnt('testCases/BNPs/PtAu20THL12_286/PtAu20THL12S2min.0.xyz', vis=False, outDir=OUTPUT_DIR, exePath=FASTBC, gridNum=1024, numPoints=300, writeBox=False)

    bcDimNPdict = {}
    boxCntDimVXPrev, boxCntDimEXPrev = 0.0, 0.0
    #for testCase in natsorted(testCases):
    for testCase in testCases:
        # if f"min.0.xyz" not in testCase: continue  # Debugging
        if 'PtAu20THL12_286' not in testCase: continue
        testCaseName = testCase.split('/')[-1].split('.')
        npName = f"{testCaseName[0]}{testCaseName[1]}_{sys.argv[2]}"
        boxCntResults = runBoxCnt(testCase, radType, radMult, calcBL, findSurfAlg, alphaMult, bulkCN, 
                                  npName, trimLen, minSample, confLvl,
                                  rmInSurf, vis, figType, saveFig, showPlot, verbose,
                                  voxelSurf, numPoints, gridNum, FASTBC, genPCD, 
                                  exactSurf, minLenMult, maxLenMult, numCPUs, numBoxLen, bufferDist, writeBox)
        r2scoreVX, boxCntDimVX, boxCntDimCIVX, minMaxLensVX, r2scoreEX, boxCntDimEX, boxCntDimCIEX, minMaxLensEX = boxCntResults
        bcDimDict = {'DBoxVX': boxCntDimVX, 'lowCIVX': boxCntDimCIVX[0], 'upCIVX': boxCntDimCIVX[1], 'R2VX': r2scoreVX, 'minLenVX': minMaxLensVX[0], 'maxLenVX': minMaxLensVX[1], 
                       'DBoxEX': boxCntDimEX, 'lowCIEX': boxCntDimCIEX[0], 'upCIEX': boxCntDimCIEX[1], 'R2EX': r2scoreEX, 'minLenEX': minMaxLensEX[0], 'maxLenEX': minMaxLensEX[1]} 
        bcDimNPdict[npName] = bcDimDict
    bcDimDF = pd.DataFrame.from_dict(bcDimNPdict, orient='index')
    bcDimDF['NPname'] = bcDimDF.index
    with open(f"BNPs{sys.argv[2]}.pickle", 'wb') as f: pickle.dump(bcDimDF, f)
    # print(bcDimDF)
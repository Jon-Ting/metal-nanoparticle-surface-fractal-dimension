from os import getcwd, listdir
from os.path import isdir
import pickle
import sys

from natsort import natsorted
import numpy as np
import pandas as pd

from sphractal.boxCnt import runBoxCnt


# Important variables
testCasePath = 'testCases/facets'
voxelSurf, exactSurf = True, False
numPoints, numBoxLen, minSample = 2500, 10, 5
radType, alphaMult = 'atomic', 2.0
minLenMult, maxLenMult = 0.15, 2.00

# General variables
trimLen, calcBL, confLvl = False, False, 95
vis, figType, rmInSurf, showPlot, writeBox, saveFig, verbose = True, 'paper', True, False, True, True, True
gridNum, findSurfAlg, genPCD = 1024, 'numNeigh', False
bufferDist, numCPUs = 5.0, int(sys.argv[1])

PROJECT_DIR = getcwd()
FASTBC = f"{PROJECT_DIR}/fastbc/3DbinImBCcpu.exe"
OUTPUT_DIR = f"{PROJECT_DIR}/outputs"


if __name__ == '__main__':

    if verbose: 
        print('Collecting test case xyz paths...')
    testCases = []
    for fName in listdir(testCasePath):
        if not isdir(f"{testCasePath}/{fName}") and fName.endswith('xyz'):
            testCases.append(f"{testCasePath}/{fName}")

    if verbose: 
        print('Running once for JIT compilation...')
    _ = runBoxCnt('testCases/BNPs/PtAu20THL12_286/PtAu20THL12S2min.0.xyz', vis=False, outDir=OUTPUT_DIR, exePath=FASTBC, findSurfAlg='numNeigh', gridNum=1024, numPoints=300, writeBox=False)

    bcDimNPdict = {}
    boxCntDimVXPrev, boxCntDimEXPrev = 0.0, 0.0
    for testCase in natsorted(testCases):
        if f"Pd31CU.xyz" not in testCase: continue  # Debugging
        npName = f"{testCase.split('/')[-1].split('.')[0]}_{sys.argv[2]}"
        print(npName)
        boxCntResults = runBoxCnt(testCase, radType, calcBL, findSurfAlg, alphaMult, 
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
    with open(f"facets{sys.argv[2]}.pickle", 'wb') as f: pickle.dump(bcDimDF, f)
    print(bcDimDF)

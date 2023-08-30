from os import getcwd, listdir
from os.path import isdir
import pickle
import sys
from time import time

from natsort import natsorted
import numpy as np
import pandas as pd

from sphractal.boxCnt import runBoxCnt


testCaseName, voxelSurf, exactSurf = 'PdNPs', False, True
PROJECT_DIR = getcwd()
FASTBC = f"{PROJECT_DIR}/sphractal/src/fastbc/3DbinImBCcpu.exe"
OUTPUT_DIR = f"{PROJECT_DIR}/outputsParamTune"

# To be tuned
# General
radTypeAlphaMult = [('atomic', 2.0), ('metallic', 2.5)]
# VX representation
gridNum = 1024  # Largest affordable resolution
numPointsAll = [2500, 5000, 10000, 20000, 40000]
numBoxLenMinSample = [(10, 3), (10, 4), (10, 5), (10, 6), (10, 7)]
# EX representation
bufferDists = [5.0, 6.0, 7.0, 8.0, 9.0]
minMaxLenMults = [(0.1, 1), (0.1, 2)]
numBoxLenMinSample = [(10, 3), (15, 5), (20, 6), (25, 8), (30, 10),  
                      #(10, 4), (15, 6), (20, 8), (25, 10), (30, 12), 
                      (10, 5), (15, 7), (20, 10), (25, 12), (30, 15)]
# Less relevant general variables
trimLen = True
vis, figType, rmInSurf, showPlot, writeBox, saveFig, verbose = True, 'paper', True, False, True, True, True
calcBL, findSurfAlg, numCPUs, confLvl, genPCD = False, 'alphaShape', int(sys.argv[1]), 95, False

# Overwrite defaults
numPoints = 10000
bufferDists = [5.0]
#minMaxLenMults = [(0.25, 1)]
#numBoxLenMinSample = [(10, 5)]


if __name__ == '__main__':

    if verbose: 
        print('Collecting test case xyz paths...')
    testCases = []
    for NPname in listdir(f"testCases/{testCaseName}"):
        testCases.append(f"testCases/{testCaseName}/{NPname}")

    if verbose: 
        print('Running once for JIT compilation...')
    _ = runBoxCnt('testCases/BNPs/PtAu20THL12_286/PtAu20THL12S2min.0.xyz', vis=False, outDir=OUTPUT_DIR, exePath=FASTBC, gridNum=512, numPoints=300, writeBox=False)

    paramTuneNPdict = {}
    for testCase in natsorted(testCases):

        if 'S1' not in testCase: continue  # Debugging
        radType, alphaMult = radTypeAlphaMult[0][0], radTypeAlphaMult[0][1]

        for bufferDist in bufferDists:
            for (minLenMult, maxLenMult) in minMaxLenMults:
                for (numBoxLen, minSample) in numBoxLenMinSample:
                    npName = f"{testCase.split('/')[-1][2:4]}_{minLenMult}_{maxLenMult}_{numBoxLen}_{minSample}_EX"
                    print('\n\n', npName)
                    start = time()
                    boxCntResults = runBoxCnt(testCase, radType, calcBL, findSurfAlg, alphaMult, 
                                              npName, trimLen, minSample, confLvl,
                                              rmInSurf, vis, figType, saveFig, showPlot, verbose,
                                              voxelSurf, numPoints, gridNum, FASTBC, genPCD, 
                                              exactSurf, minLenMult, maxLenMult, numCPUs, numBoxLen, bufferDist, writeBox)
                    # Recommend maxLenMult = 2 for 'full'
                    duration = time() - start
                    r2scoreVX, boxCntDimVX, boxCntDimCIVX, minMaxLensVX, r2scoreEX, boxCntDimEX, boxCntDimCIEX, minMaxLensEX = boxCntResults
    
                    paramTuneDict = {'DBoxVX': boxCntDimVX, 'lowCIVX': boxCntDimCIVX[0], 'upCIVX': boxCntDimCIVX[1], 'R2VX': r2scoreVX, 'minLenVX': minMaxLensVX[0], 'maxLenVX': minMaxLensVX[1], 
                                   'DBoxEX': boxCntDimEX, 'lowCIEX': boxCntDimCIEX[0], 'upCIEX': boxCntDimCIEX[1], 'R2EX': r2scoreEX, 'minLenEX': minMaxLensEX[0], 'maxLenEX': minMaxLensEX[1], 
                                   'radType': radType, 'findSurfAlg': findSurfAlg, 'alphaMult': alphaMult, 'trimLen': trimLen, 'minSample': minSample, 'confLvl': confLvl, 'rmInSurf': rmInSurf,
                                   'numPoints': numPoints, 'gridNum': gridNum, 
                                   'minLenMult': minLenMult, 'maxLenMult': maxLenMult, 'numBoxLen': numBoxLen, 'bufferDist': bufferDist, 'duration': duration}
                    paramTuneNPdict[npName] = paramTuneDict
                    print(f"DBoxVX: {boxCntDimVX:.4f}, R2VX: {r2scoreVX:.4f}, DBoxEX: {boxCntDimEX:.4f}, R2EX: {r2scoreEX:.4f}, duration: {duration:.4f}")

    paramTuneDF = pd.DataFrame.from_dict(paramTuneNPdict, orient='index')
    with open(f"paramTune{sys.argv[2]}.pickle", 'wb') as f: pickle.dump(paramTuneDF, f)
    print(paramTuneDF)

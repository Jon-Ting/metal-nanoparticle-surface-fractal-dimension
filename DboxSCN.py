from os import getcwd, listdir
from os.path import isdir
import pickle
import sys

from natsort import natsorted
import numpy as np
import pandas as pd

from sphractal.boxCnt import runBoxCnt


testCasePath = '/mnt/c/Users/ASUS/Documents/PhD/Workstation/Data/Characterised/Palladium_Nanoparticle_Data_Set'
voxelSurf, exactSurf = True, True
numPoints, numBoxLen, minSample = 300, 10, 5
radType, alphaMult = 'atomic', 2.0
minLenMult, maxLenMult = 0.25, 1

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
    for fName in listdir(testCasePath):
        if not isdir(f"{testCasesPath}/{fName}") and fName.endswith('xyz'):
            testCases.append(f"testCases/{testCaseName}/{NPname}")

    if verbose: 
        print('Running once for JIT compilation...')
    _ = runBoxCnt('testCases/BNPs/PtAu20THL12_286/PtAu20THL12S2min.0.xyz', vis=False, outDir=OUTPUT_DIR, exePath=FASTBC, gridNum=1024, numPoints=300, writeBox=False)

    fracDimNPdict = {}
    for testCase in natsorted(testCases):
        # if f"Pd{sys.argv[2]}" not in testCase: continue  # Debugging
        boxCntResults = runBoxCnt(testCase, radType, calcBL, findSurfAlg, alphaMult, 
                          npName, trimLen, minSample, confLvl,
                          rmInSurf, vis, figType, saveFig, showPlot, verbose,
                          voxelSurf, numPoints, gridNum, FASTBC, genPCD, 
                          exactSurf, minLenMult, maxLenMult, numCPUs, numBoxLen, bufferDist, writeBox)
        # Recommend maxLenMult = 2 for 'full'
        r2scoreVX, boxCntDimVX, boxCntDimCIVX, minMaxLensVX, r2scoreEX, boxCntDimEX, boxCntDimCIEX, minMaxLensEX = boxCntResults
        fracDimDict = {'DBoxVX': boxCntDimVX, 'lowCIVX': boxCntDimCIVX[0], 'upCIVX': boxCntDimCIVX[1], 'R2VX': r2scoreVX, 'minLenVX': minMaxLensVX[0], 'maxLenVX': minMaxLensVX[1], 
                       'DBoxEX': boxCntDimEX, 'lowCIEX': boxCntDimCIEX[0], 'upCIEX': boxCntDimCIEX[1], 'R2EX': r2scoreEX, 'minLenEX': minMaxLensEX[0], 'maxLenEX': minMaxLensEX[1]} 
        fracDimNPdict[testCase] = fracDimDict
    fracDimDF = pd.DataFrame.from_dict(fracDimNPdict, orient='index')
    fracDimDF['NPname'] = fracDimDF.index
    # fracDimDF['NPshape'] = fracDimDF['NPname'].apply(lambda x: x[-8:-6])
    with open(f"DboxSCN.pickle", 'wb') as f: pickle.dump(fracDimDF, f)
    # with open(f"DBox.pickle", 'rb') as f: fracDimDF = pickle.load(f)
    print(fracDimDF)

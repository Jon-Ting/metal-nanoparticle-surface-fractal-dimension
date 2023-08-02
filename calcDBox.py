from os import listdir
from os.path import isdir

from natsort import natsorted
import numpy as np
import pandas as pd

from sphractal import runBoxCnt


testCaseName, boxLenRange, runPointCloudBoxCnt, runExactSurfBoxCnt = 'PtAu20THL12_286', 'Trimmed', True, False  # CoPd150CUL10_307063, PtAu20THL12_286
vis, writeBox, rmInSurf, verbose, calcAvgDuration, findSurfOption = True, True, True, True, False, 'alphaShape'
boxLenConc, atomConc, boxConc = False, True, False
GRID_NUM = 1024  # 4056 sufficient for < 3500 points for disPd1, Max 1048576 for laptop before RAM runs out
NUM_SPHERE_POINT = 300  # Number of points to be fitted onto each atomic sphere
ALPHA_MULT = 2.5
PROJECT_DIR = '/mnt/c/Users/ASUS/Documents/PhD/Workstation/PaperDrafts/NP_Frac_Dim/DOCS'
BIN_IM_BC_EXE_DIR = f"{PROJECT_DIR}/sphractal/bin"
OUTPUT_DIR = f"{PROJECT_DIR}/testOutputs"


@utils.estDuration
def calcBoxCntDim(atoms, maxDimDiff, eleSet, atomPosMinMax, testCase,
                  findSurfOption='alphaShape', calcBL=False):
    """For running on NCI HPC facilities (Gadi) where data are stored."""
    minAtomRad, maxAtomRad, minRadEle, maxRadEle = 10.0, 0.0, '', ''
    for atomEle in eleSet:
        if constants.ATOMIC_RAD_DICT[atomEle] < minAtomRad:
            minAtomRad, minRadEle = constants.ATOMIC_RAD_DICT[atomEle], atomEle
        if constants.ATOMIC_RAD_DICT[atomEle] > maxAtomRad:
            maxAtomRad, maxRadEle = constants.ATOMIC_RAD_DICT[atomEle], atomEle
    utils.findNN(atoms, atomPosMinMax, maxAtomRad, calcBL=calcBL)
    utils.findSurf(atoms, option=findSurfOption, alpha=2*constants.ATOMIC_RAD_DICT[minRadEle])
    minBoxLen, maxBoxLen = minAtomRad / 4, minAtomRad
    # TODO: Remove @estDuration first!
    scaleChange, countChange = boxCnt.getSphereBoxCnts(atomList, maxDimDiff, (minBoxLen, maxBoxLen),
                                                       atomPosMinMax[:3], OUTPUT_DIR, testCase,
                                                       rmInSurf=True, writeBox=True, verbose=False,
                                                       boxLenConc=False, atomConc=True,
                                                       boxConc=False)
    r2score, boxCntDim, slopeCI = boxCnt.findSlope(scaleChange, countChange,
                                                   OUTPUT_DIR, testCase, boxLenRange='Trimmed',
                                                   visReg=True, saveFig=True, showPlot=False)
    avgBondLens = np.array([atom.avgBondLen for atom in atoms]) if calcBL else None
    return r2score, boxCntDim, slopeCI, avgBondLens


if __name__ == '__main__':
    testCases = []
    for NPname in listdir(f"testCases/{testCaseName}"):
        if not isdir(f"testCases/{testCaseName}/{NPname}"): testCases.append(f"testCases/{testCaseName}/{NPname}")
        else: testCases.extend(natsorted([f"testCases/{NPname}/{i}" for i in listdir(f"testCases/{NPname}")]))

    # For first just-in-time compilation
    _ = runBoxCnt('testCases/PtAu20THL12_286/PtAu20THL12S2min.0.xyz', runPointCloudBoxCnt=False, vis=False, writeBox=False)

    fracDimNPdict = {}
    boxCntDimPrev = 0.0
    for testCase in testCases:
        # if '000.xyz' not in testCase: continue  # Debugging
        boxCntResults = runBoxCnt(testCase, radType, calcBL, findSurfOption, alphaMult, 
                                  OUTPUT_DIR, lenRange, minSampleNum, confLvl,
                                  rmInSurf, vis, figType, saveFig, showPlot, verbose,
                                  runPointCloudBoxCnt, numSpherePoint, gridNum, FASTBC_EXE, genPCD, 
                                  runExactSphereBoxCnt, 0.25, 1, numBoxLenSample, writeBox)
        boxCntDimPC, boxCntDimCIPC, r2scorePC, boxCntDimES, boxCntDimCIES, r2scoreES = boxCntResults
        print(f"{testCase}\t\tD_Box: {boxCntDim:.4f} [{boxCntDimCI[0]:.4f}, {boxCntDimCI[1]:.4f}]\t\tR2: {r2score:.4f}"
              f"\t\tTime: {duration:.4f}\t\tDiff: {abs(boxCntDim - boxCntDimPrev):.4f}")
        boxCntDimPrev = boxCntDim
        fracDimDict = {'boxCntDim': boxCntDim, 'lowCI': boxCntDimCI[0], 'upCI': boxCntDimCI[1], 'R2': r2score,
                       't': duration}
        fracDimNPdict[testCase] = fracDimDict

    fracDimDF = pd.DataFrame.from_dict(fracDimNPdict, orient='index')
    fracDimDF['NPname'] = fracDimDF.index
    fracDimDF['NPshape'] = fracDimDF['NPname'].apply(lambda x: x[:2])
    fracDimDF['NPtemp'] = fracDimDF['NPname'].apply(lambda x: x[-3:])
    # TH = 4 {111}
    # CU = 6 {100}
    # OT = 8 {111}
    # RD = 12 {110}
    # TO = 8 {111}, 6 {100}
    # CO = 8 {111}, 6 {100}
    # DH = 10 {111}, 5 {100}  # Ino
    # IC = 20 {111}
    numFaceOrderRank = {'TH': 1, 'CU': 2, 'OT': 3, 'RD': 4, 'TO': 5, 'CO': 6, 'DH': 7, 'IC': 8, 'di': 9}
    fracDimDF['orderRank'] = fracDimDF['NPshape'].map(numFaceOrderRank)
    fracDimDF.sort_values(by='orderRank')
    fracDimIdealNPDF = fracDimDF[fracDimDF['NPtemp'] == '000']
    fracDimIdealNPDF.groupby(by='NPshape').mean().sort_values(by='orderRank')
    fracDimHeatedNPDF = fracDimDF[fracDimDF['NPtemp'] == '323']
    fracDimHeatedNPDF.groupby(by='NPshape').mean().sort_values(by='orderRank')
    
    import pickle
    with open(f"{boxLenRange}.pickle", 'wb') as f: pickle.dump(fracDimDF, f)
    with open(f"{boxLenRange}.pickle", 'rb') as f: fracDimDF = pickle.load(f)
    print(fracDimDF)
    
    if verbose:
        print(f"Change in\n\tScale: {scaleChange}\n\tCounts: {countChange}")
        print(f"Coefficient of determination (R2): {r2score:.3f}\
              \nEstimated box-counting dimension (D_Box): {boxCntDim:.3f}\
              \n{CONF_INT_PERC}% Confidence interval: [{boxCntDimCI[0]:.3f}, {boxCntDimCI[1]:.3f}]")

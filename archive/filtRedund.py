"""
Help on package filtRedund:

NAME
    filtRedund

DESCRIPTION
    Module to filter redundant snapshots from molecular dynamics simulations of metal nanoparticles
    ===============================================================================================
    Implemented:
    - Euclidean distance RMSD
    - Radial Distribution Function Mutual Information (didn't take element type into account, could use NCPac as well)
    - Q6Q6 average distribution (Requires NCPac to be run first, which defeats the purpose of filtering redundant frames prior to feature extraction, thus omitted)
    - FD values
"""

import os
import warnings

from natsort import natsorted
from pyinform.mutualinfo import mutual_info
from rdfpy import rdf
from scipy.stats import ks_2samp, cramervonmises_2samp, anderson_ksamp, wasserstein_distance, gaussian_kde
from scipy.special import rel_entr
from sklearn.feature_selection import mutual_info_regression

from FDestimation import *  # Includes numpy and time


# Control printout to ease debugging, to be commented out eventually
warnings.filterwarnings(action='ignore', category=RuntimeWarning)  # Addresses warnings from 2-sample Kolmogorov-Smirnov test regarding switching to asymp from exact for p-value calculation
warnings.filterwarnings(action='ignore', category=UserWarning)  # Addresses warnings from k-sample Anderson-Darling test regarding flooring/capping p-value

# eleComb = 'AuPt'
# NP_DIRS_PATH = f"/scratch/q27/jt5911/SimAnneal/{eleComb}"
NP_DIRS_PATH = "./testCases/AuPt40CORCS_3871"
PVAL_THRESH = 0.05
RDF_RCUTOFF_MULTIPLIER = 0.9  # Minimise noise at edge particles
RDF_NUM_SAMPLES = 1000  # Balance between computational cost and statistical test power

# Fine-tuned based on AuPd20CUL10, AuPd30ICRCCS_923, CoPd40TOL12_4033
EUC_DIST_THRESH = 0.05  # 0.003-0.063
BOX_CNT_DIM_DIFF_THRESH = 0.01  # 0.0061-0.0229
WS_THRESH = 0.0001  # 0.00004-0.00020
HL_THRESH = 0.05  # 0.03-0.09
KL_THRESH = 0.2  # 0.16-0.27
JS_THRESH = 0.005  # 0.002-0.008
MI_THRESH = 1.0  # 0.86-1.29
# RDF_DIFF_THRESH = 0.5  # 0.12 (visually undetectable) to 0.6 (1 atom moved)
# BC_THRESH = 0.89
# Q6Q6_BOND_DIFF_THRESH = 0.001


def alignNP(P, Q):
    '''
    Align P with Q (reference structure) using Kabsch algorithm
    Algorithm explanation: https://en.wikipedia.org/wiki/Kabsch_algorithm
    Reference code: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
    '''
    # Translation of centroids of both structures to system origin
    P_cent, Q_cent = P - P.mean(axis=0), Q - Q.mean(axis=0)
    C = np.dot(P_cent.T, Q_cent)  # Covariance matrix
    V, S, W = np.linalg.svd(C)  # Singular value decomposition
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0  # Proper/improper rotation
    if d: S[-1], V[:, -1], = -S[-1], -V[:, -1] # Antialigns of the last singular vector
    U = np.dot(V, W)  # Rotation matrix
    P_align = np.dot(P_cent, U) + P.mean(axis=0)
    return P_align
    

@estDuration
def calcEucRMSD(coordList1, coordList2):
    eucRMSD = np.sqrt(np.mean(np.square(np.array(coordList1) - np.array(coordList2))))
    return eucRMSD


def getCDF(p, covFactor=0.1):
    cdf = gaussian_kde(p)
    cdf.covariance_factor = lambda: covFactor
    cdf._compute_covariance()
    return cdf


def calcBCdist(p, q):
    # https://github.com/EricPWilliamson/bhattacharyya-distance/blob/master/bhatta_dist.py
    pCDF, qCDF = getCDF(p), getCDF(q)
    minVal, maxVal, numSteps = min(np.concatenate((p, q))), max(np.concatenate((p, q))), 100
    xs = np.linspace(minVal, maxVal, numSteps)
    return sum([np.sqrt(pCDF(x)[0]*qCDF(x)[0])*(maxVal-minVal)/numSteps for x in xs])


def calcHLdist(p, q):
    return np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)


def calcKLdiv(p, q):
    return np.sum(rel_entr(p+10e-128, q+10e-128))  # rel_entr == p*np.log(p/q)


def calcMI(x, y):
    # To calculate MI, joint distribution (pair distribution) is needed, not just marginal distributions
    return mutual_info(x, y)


@estDuration
def calcRDFdiff(coordList1, coordList2, minMaxDimDiff, showPlot=False):
    rdfInt = round(minMaxDimDiff/2 * RDF_RCUTOFF_MULTIPLIER / RDF_NUM_SAMPLES, 3)
    gR1, radii1 = rdf(coordList1, dr=rdfInt)
    gR2, radii2 = rdf(coordList2, dr=rdfInt)
    while np.isnan(gR1).any() or np.isnan(gR2).any():
        gR1, gR2 = gR1[:-1], gR2[:-1]
    A = anderson_ksamp([gR1, gR2], midrank=True)  # H0: F=G
    D = ks_2samp(gR1, gR2, alternative='two-sided')  # 0 (same) < D < 1 (different)
    W = cramervonmises_2samp(gR1, gR2, method='asymptotic')  # 'exact' option suffers from combinatorial explosion
    # if showPlot:
    #     import matplotlib.pyplot as plt
    #     plt.plot(radii1, gR1);
    #     plt.plot(radii2, gR2);
    #     plt.xlabel('r');
    #     plt.ylabel('g(r)');
    #     plt.legend(['NP1', 'NP2']);
    #     plt.show();  # plt.savefig('test.png');
    if len(gR1) != len(gR2):  # Cut the last points from the slightly longer g_r
        if len(gR1) < len(gR2): gR2 = gR2[:len(gR1)]
        else: gR1 = gR1[:len(gR2)]
    gR1norm, gR2norm = gR1/gR1.sum(), gR2/gR2.sum()
    bcDist = calcBCdist(gR1norm, gR2norm) # Bhattacharyya distance
    wsDist = wasserstein_distance(gR1norm, gR2norm, u_weights=None, v_weights=None)  # Wasserstein distance
    hlDist = calcHLdist(gR1norm, gR2norm)  # Hellinger distance
    klDiv = calcKLdiv(gR1norm, gR2norm)  # Kullback-Leibler divergence, asymmetrical
    m = 0.5*(gR1norm+gR2norm)
    jsDiv = 0.5*calcKLdiv(gR1norm, m) + 0.5*calcKLdiv(gR2norm, m)  # Jensen-Shannon divergence, symmetrical
    miScore = calcMI(gR1, gR2)  # Mutual information
    return A.statistic, A.pvalue, D.statistic, D.pvalue, W.statistic, W.pvalue, hlDist, klDiv, jsDiv, miScore, bcDist, wsDist  #, np.sqrt(np.sum(np.square(gR1 - gR2)))


# def readQ6Q6(npName):
#     # Need to run NCPac.exe first, might not be feasible for large nanoparticles
#     q6q6bondList = []
#     with open(f"{NP_DIRS_PATH}/{npName}/ov_Q6Q6.xyz", 'r') as f:
#         for (i, line) in enumerate(f):
#             if i < 2: continue
#             q6q6bondNum = float(line.split()[-1])
#             q6q6bondList.append(q6q6bondNum)
#     return q6q6bondList


# def calcQ6Q6bondDiff(npName1, npName2):
#     q6q6bonds1, q6q6bonds2 = readQ6Q6(npName1), readQ6Q6(npName2) 
#     q6q6Diff = np.sqrt(np.sum(np.square(np.array(q6q6bonds1) - np.array(q6q6bonds2))))
#     return q6q6Diff


def runFilter(calcBL=False):
    foundRedund = False
    NPfiltIdxs = set()
    sortedDirs = natsorted(os.listdir(NP_DIRS_PATH))
    print(f"Filtering redundant frames from {sortedDirs[0].split('min')[0]} nanoparticle simulation...")
    eucTs, rdfTs, bcdTs, blcnTs = [], [], [], []  # q6q6Ts = []
    for (i, npName1) in enumerate(sortedDirs):
        # Compute box-counting dimension for 1st NP
        if i not in NPfiltIdxs:
            atomList1, maxDimDiff1, eleSet1, atomPosMinMax1 = readXYZ(f"{NP_DIRS_PATH}/{npName1}")
            coordList1 = np.array([[atom.X, atom.Y, atom.Z] for atom in atomList1])
            (r2dim1, boxCntDim1, CIdim1, avgBondLen1, coordNum1), bcdT1 = calcBoxCntDim(atomList1, maxDimDiff1, eleSet1, atomPosMinMax1, calcBL=calcBL)
            (avgBondLen1, coordNum1), blcnT1 = calcBLCN(atomList1, eleSet1, atomPosMinMax1)

        # Determine the second NP to compare to
        j = j + 1 if foundRedund else i + 1
        if j >= len(sortedDirs): break
        npName2 = sortedDirs[j]
        # if j != i+1: continue  # j!=i+1 (consecutive) / j>=i (all other)
        print(f"  Comparing {npName1} with {npName2}...")

        # Compute box-counting dimension for 2nd NP
        atomList2, maxDimDiff2, eleSet2, atomPosMinMax2 = readXYZ(f"{NP_DIRS_PATH}/{npName2}")
        coordList2 = np.array([[atom.X, atom.Y, atom.Z] for atom in atomList2])
        coordList2 = alignNP(coordList2, coordList1)  # Doesn't impact FD estimation
        (r2dim2, boxCntDim2, CIdim2, avgBondLen2, coordNum2), bcdT2 = calcBoxCntDim(atomList2, maxDimDiff2, eleSet2, atomPosMinMax2, calcBL=calcBL)
        (avgBondLen2, coordNum2), blcnT2 = calcBLCN(atomList2, eleSet2, atomPosMinMax2)

        # Compute and report all difference metrics
        eucRMSD, eucT = calcEucRMSD(coordList1, coordList2)
        (rdfStatA, rdfPValA, rdfStatD, rdfPValD, rdfStatW, rdfPValW, rdfHL, rdfKL, rdfJS, rdfMI, rdfBC, rdfWS), rdfT = calcRDFdiff(coordList1, coordList2, min(maxDimDiff1, maxDimDiff2))
        # q6q6Diff = calcQ6Q6bondDiff(npName1, npName2)
        boxCntDimDiff = abs(boxCntDim2 - boxCntDim1)
        avgBondLenMI = mutual_info_regression(avgBondLen1.reshape((-1, 1)), avgBondLen2, random_state=777)[0]
        coordNumMI = mutual_info_regression(coordNum1.reshape((-1, 1)), coordNum2, random_state=777)[0]

        eucTs.append(eucT)
        rdfTs.append(rdfT)
        bcdTs.append(bcdT1 + bcdT2)
        blcnTs.append(blcnT1 + blcnT2)
        print(f"    eucRMSD: {eucRMSD:.3f}")
        # print(f"    q6q6Diff: {q6q6Diff}")
        print(f"    boxCntDimDiff: |{boxCntDim1:.4f} - {boxCntDim2:.4f}| = {boxCntDimDiff:.4f}")
        print(f"    RDF 2-sample Anderson-Darling test statistic: {rdfStatA:.4f}; p-value: {rdfPValA:.3f}")
        print(f"    RDF 2-sample Kolmogorov-Smirnov test statistic: {rdfStatD:.4f}; p-value: {rdfPValD:.3f}")
        print(f"    RDF 2-sample Cramer-von Mises test statistic: {rdfStatW:.4f}; p-value: {rdfPValW:.3f}")
        print(f"    RDF Bhattacharyya distance : {rdfBC:.2f}")  # TODO: switch off
        print(f"    RDF Wasserstein distance : {rdfWS:.7f}")
        print(f"    RDF Hellinger distance : {rdfHL:.3f}")
        print(f"    RDF Kullback-Leibler divergence : {rdfKL:.2f}")
        print(f"    RDF Jensen-Shannon divergence : {rdfJS:.5f}")
        print(f"    RDF Mutual information : {rdfMI:.2f}")
        print(f"    avgBondLen Mutual information : {avgBondLenMI:.2f}")
        print(f"    coordNum Mutual information : {coordNumMI:.2f}")  # TODO: switch off

        # Determine if 1st NP is redundant
        if eucRMSD < EUC_DIST_THRESH and boxCntDimDiff < BOX_CNT_DIM_DIFF_THRESH and rdfPValA > PVAL_THRESH and rdfPValD > PVAL_THRESH and rdfPValW > PVAL_THRESH and rdfWS < WS_THRESH and rdfHL < HL_THRESH and rdfKL < KL_THRESH and rdfJS < JS_THRESH and rdfMI > MI_THRESH:
            print(f"    {npName1} redundant!")
            NPfiltIdxs.add(i)
            foundRedund = True
            continue
        foundRedund = False
    print(f"Average duration (s):\n  RMSD {sum(eucTs)/len(eucTs):.6f}\n  RDF {sum(rdfTs)/len(rdfTs):.6f}\n  BCDim {sum(bcdTs)/len(bcdTs):.6f}\n  BondLenCoordNum {sum(blcnTs)/len(blcnTs):.6f}")
    return NPfiltIdxs


if __name__ == '__main__':
    NPfiltIdxs = runFilter(calcBL=False)
    print(NPfiltIdxs)

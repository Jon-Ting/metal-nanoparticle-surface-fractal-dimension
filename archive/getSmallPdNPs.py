from os import listdir
from shutil import copy


smallNPnames = []
for fName in listdir('.'):
    if 'xyz' in fName:
        with open(fName, 'r') as f1:
            numAtom = int(f1.readline())
            if numAtom <= 1000:
                print(fName)
                smallNPnames.append(fName.split('.')[0])
                copy(fName, './smallPdNPs/')

with open('Pd_nanoparticle_dataset.csv', 'r') as f2, open('smallPdNPs.csv', 'w') as f3:
    for (i, line) in enumerate(f2):
        if i == 0:
            continue
        if line.split(',')[0].zfill(6) in smallNPnames:
            print(i)
            f3.write(line)

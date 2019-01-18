import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import csv
import pandas as pd
from collections import defaultdict
import pysam
import os
import pysamstats
from collections import Counter


def load_rdata(rdata):
    pandas2ri.activate()
    robjects.r.load(rdata)

    return robjects.r


def load_bininfo(bininfoData):
    with open(bininfoData, 'r') as infile:
        csvList = csv.reader(infile)
        header = list()
        elems = list()
        for i, line in enumerate(csvList):
            if not i == 0:
                elems.append(line)
            else:
                header = line

        bininfo = pd.DataFrame(elems, columns=header)

    bininfo.rename(columns={'': "binName"}, inplace=True)
    bininfo['binorder'] = range(1, 61928)

    return bininfo


def load_sam(file):
    bincount = defaultdict(int)
    exceptions = ['*', 'chrM']
    with open(file, 'r') as infile:
        for line in infile:
            elem = line.strip().split("\t")  # could be changed
            if elem[2] not in exceptions:
                chrPos = elem[2] + "_" + str( (int(elem[3]) -1) // 50000)
                bincount[chrPos] += 1

    return bincount


def load_counts(file):
    bincount = {}
    with open(file, 'r') as infile:
        for line in infile:
            line = line.strip()
            elem = line.split("\t")  # could be changed

            chrPos = elem[0]
            count = int(elem[1])

            bincount[chrPos] = count

    return bincount


def load_bam(file):
    print(file)
    bam = pysam.AlignmentFile(file)
    homo_sapiens = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                    'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                    'chrX', 'chrY']

    bincount = defaultdict(int)

    for i, read in enumerate(bam):
        chrom = read.reference_name
        #if chrom not in exceptions:
        chrPos = str(chrom) + "_" + str(int(read.pos) // 50000)   # 0-based
        bincount[chrPos] += 1

    exceptions1 = [chrom for chrom in bincount.keys() if chrom.split("_")[0] not in homo_sapiens]
    exceptions2 = [chrom for chrom in bincount.keys() if len(chrom.split("_")) != 2]
    exceptions3 = set(exceptions1).union(set(exceptions2))

    for exception in exceptions3:
        bincount.pop(exception)

    print("bam file completely loaded : {}".format(os.path.split(file)[1]))
    return bincount
import csv
import numpy as np
import pandas as pd
import random
import glob
import os
import time
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt

"""
1. load files
2. calculate chrY counts / all counts 
"""


def fraction_chrY(path):
    fraction_list = []
    filename_list = []

    random.seed(time.time())

    sampling=random.sample(glob.glob(path), 10)

    for file in sampling:
        allcounts = 0
        chrYcounts = 0
        filename_list.append(os.path.basename(file))

        with open(file, 'r') as infile:
            for line in infile:
                line = line.strip()
                elem = line.split("\t")
                allcounts += int(elem[1])
                if 'chrY' in elem[0]: chrYcounts += int(elem[1])

        fraction_list.append(chrYcounts / allcounts)

    return fraction_list, filename_list


def percentage_chrY(path):
    fraction_list, _ = fraction_chrY(path)
    perc_chrY = np.mean(fraction_list, dtype=np.float64)

    return perc_chrY


def ff_chrY(female_perc_chrY, male_perc_chrY, newtemp_test):

    ff_list = []
    test_fraction, filename_list = fraction_chrY(newtemp_test)

    #print(test_fraction)
    for test_perc_chrY in test_fraction:
        ff = (test_perc_chrY - female_perc_chrY) / (male_perc_chrY - female_perc_chrY)
        ff_list.append(ff)

    return ff_list, filename_list


if __name__ == "__main__":
    # female_1000 = percentage_chrY("/home/hyunbin/seqff/output/female/newtemp/*.newtemp")
    # male_6 = percentage_chrY("/home/hyunbin/seqff/output/adult_male/newtemp/*.newtemp")

    female_1000 = 0.0023164119684942587
    male_6 = 0.00714710817264482

    ff_list, filename_list = ff_chrY(female_1000, male_6, "/home/hyunbin/seqff/output/male/newtemp/*.newtemp")

    ff_df = pd.DataFrame(data={'Dataname': filename_list, 'FetalFraction': ff_list})
    #save result as a csv file
    #ff_df.to_csv("/home/hyunbin/ff_chrY.csv", index=False, header=False)

    # plt.hist(ff_list)
    # plt.show()
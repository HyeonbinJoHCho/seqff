import csv
import sys
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


def fraction_chrY(path, gender="female"):
    fraction_list = []
    filename_list = []

    random.seed(time.time())

    if gender == "female":
        #sampling = random.sample(glob.glob(path), 10)
        sampling = glob.glob(path)
    elif gender in ("adult_male", "male_fetus"):
        sampling = glob.glob(path)

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


#TEMPORARY
def print_chrY_percentages(fraction_list, filename_list, gender):
    tmp_df = pd.DataFrame(data= {'filename':filename_list, 'chrYfraction':fraction_list})
    if gender == 'female':
        tmp_df.to_csv("/home/hyunbin/female_1000_chrY_fraction.csv", index=False, header=False)
    elif gender == "adult_male":
        tmp_df.to_csv("/home/hyunbin/male_6_chrY_fraction.csv", index=False, header=False)
    elif gender == "male_fetus":
        tmp_df.to_csv("/home/hyunbin/male_fetus_1000_chrY_fraction.csv", index=False, header=False)


def percentage_chrY(path, gender):
    fraction_list, filename_list = fraction_chrY(path, gender)
    perc_chrY = np.mean(fraction_list, dtype=np.float64)

    #print_chrY_percentages(fraction_list, filename_list, gender=gender)

    return perc_chrY


def ff_chrY(female_perc_chrY, male_perc_chrY, newtemp_test):

    ff_list = []
    test_fraction, filename_list = fraction_chrY(newtemp_test)

    for test_perc_chrY in test_fraction:
        ff = (test_perc_chrY - female_perc_chrY) / (male_perc_chrY - female_perc_chrY)
        ff_list.append(ff)

    return ff_list, filename_list


if __name__ == "__main__":
    #female_1000 = percentage_chrY("/home/hyunbin/seqff/output/female/newtemp/*.newtemp", "female")
    #male_6 = percentage_chrY("/home/hyunbin/seqff/output/adult_male/newtemp/*.newtemp", "adult_male")

    female_1000 = 0.0023164119684942587
    male_6 = 0.00714710817264482
    print("female_6", female_1000)
    print("male_6", male_6)

    ff_list, filename_list = ff_chrY(female_1000, male_6, "/home/hyunbin/seqff/output/male_20190114/newtemp/*.newtemp")

    ff_df = pd.DataFrame(data={'Dataname': filename_list, 'FetalFraction': ff_list})
    #save result as a csv file
    #ff_df.to_csv("/home/hyunbin/seqff/output/male_20190114/ff_chrY.csv", index=False, header=False)

    # plt.hist(ff_list)
    # plt.show()

    print("completed")
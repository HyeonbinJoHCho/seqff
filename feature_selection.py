import pandas as pd
import numpy as np
import os, time, sys, glob
import sklearn.feature_selection as fs
import math
import traceback

from sklearn.linear_model import MultiTaskElasticNetCV, MultiTaskElasticNet

import matplotlib.pyplot as plt


class FeatureSelection():
    def __init__(self):
        ##location info
        self.female_path = glob.glob("/home/hyunbin/seqff/output/female_notfiltered/normalizedbincount/*")
        self.male_path = glob.glob("/home/hyunbin/seqff/output/male_notfiltered/normalizedbincount/*")
        self.supple_path = "/home/hyunbin/seqff/supplementary-table2.csv"
        self.alluseablebins_path = "/home/hyunbin/seqff/output/male_notfiltered/alluseablebins"

    def pearson_corr(self, all_auto_df):
        all_auto_array = np.array([all_auto_df[i].values for i in all_auto_df.columns])
        pearson_result = np.corrcoef(all_auto_array[:, 1000:])
        print("pearson result is completely calculated")

    def selection_with_corr(self):
        raise NotImplementedError

    def selection_with_elastic_net_(self, all_auto_df, all_y_df):
        start = time.time()

        all_auto_array = np.array([all_auto_df[i].values for i in all_auto_df.columns])
        all_y_array = np.array([all_y_df[i].values for i in all_y_df.columns])

        l1_ratio = 0.01  # could be changed

        regr = MultiTaskElasticNetCV(l1_ratio=l1_ratio, n_alphas=100,
                                     cv=10)
        regr.fit(all_auto_array[:, 1000:].transpose(), all_y_array[:, 1000:].transpose())
        end = time.time()

        elapsed = end - start

        print("alpha:", regr.alpha_)
        print("l1_ratio_:", regr.l1_ratio_)
        print(regr.coef_.shape)

        print("%d" % elapsed)
        return regr

    def result(self, select):
        supple = pd.read_csv(self.supple_path)
        alluseablebins = pd.read_csv(self.alluseablebins_path, header=None)

        supple['alluseablebins'] = alluseablebins

        female_df_list = []
        male_df_list = []
        for file in self.female_path:
            with open(file, 'r') as infile:
                df = pd.read_csv(infile, sep="\t", header=None)
                df.columns = ['normalizedbincount']
                df['Unnamed: 0'] = np.array(supple[supple['alluseablebins'] == True]['Unnamed: 0'])
                df['CHR'] = np.array(supple[supple['alluseablebins'] == True]['CHR'])
                df.fillna(0, inplace=True)
                female_df_list.append(df)

        for file in self.male_path:
            with open(file, 'r') as infile:
                df = pd.read_csv(infile, sep="\t", header=None)
                df.columns = ['normalizedbincount']
                df['Unnamed: 0'] = np.array(supple[supple['alluseablebins'] == True]['Unnamed: 0'])
                df['CHR'] = np.array(supple[supple['alluseablebins'] == True]['CHR'])
                df.fillna(0, inplace=True)
                male_df_list.append(df)

        female_auto_df_list = []
        male_auto_df_list = []

        for index, df in enumerate(female_df_list):
            temp_df = df[(df['CHR'] != "chrY")]
            temp_df = temp_df[(temp_df['CHR'] != "chrX")]
            female_auto_df_list.append(temp_df)

        for index, df in enumerate(male_df_list):
            temp_df = df[(df['CHR'] != "chrY")]
            temp_df = temp_df[(temp_df['CHR'] != "chrX")]
            male_auto_df_list.append(temp_df)

        female_y_df_list = []
        male_y_df_list = []

        for index, df in enumerate(female_df_list):
            # female_auto_df_list[index] = female_auto_df_list[index][(female_auto_df_list[index][0] != 0)]
            temp_df = df[(df['CHR'] == "chrY")]
            female_y_df_list.append(temp_df)
        for index, df in enumerate(male_df_list):
            # male_auto_df_list[index] = male_auto_df_list[index][(male_auto_df_list[index][0] != 0)]
            temp_df = df[(df['CHR'] == "chrY")]
            male_y_df_list.append(temp_df)

        female_auto_concat = pd.concat([df['normalizedbincount'] for df in female_auto_df_list], axis=1)
        male_auto_concat = pd.concat([df['normalizedbincount'] for df in male_auto_df_list], axis=1)

        all_auto_df = pd.concat([female_auto_concat, male_auto_concat], axis=1)
        all_auto_df.index = male_auto_df_list[0]['Unnamed: 0']
        all_auto_df.columns = ["sample{}".format(i) for i in range(1, len(all_auto_df.columns) + 1)]
        all_auto_df = all_auto_df.transpose()

        female_y_concat = pd.concat([df['normalizedbincount'] for df in female_y_df_list], axis=1)
        male_y_concat = pd.concat([df['normalizedbincount'] for df in male_y_df_list], axis=1)

        all_y_df = pd.concat([female_y_concat, male_y_concat], axis=1)
        all_y_df.index = male_y_df_list[0]['Unnamed: 0']
        all_y_df.columns = ["sample{}".format(i) for i in range(1, len(all_y_df.columns) + 1)]
        all_y_df = all_y_df.transpose()

        if select == 'pearson':
            corr = self.pearson_corr(all_auto_df)
            #TODO
            #save result to csvfile

        elif select == 'enet':
            enet_result = self.selection_with_elastic_net_(all_auto_df, all_y_df)

        elif select == '':
            raise NotImplementedError



if __name__ == "__main__":
    featureSelection = FeatureSelection()
    featureSelection.result(select='enet')

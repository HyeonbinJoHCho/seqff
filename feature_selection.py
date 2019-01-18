import pandas as pd
import numpy as np
import os, time, sys, glob
import sklearn.feature_selection as fs
import math
import traceback
import csv
import itertools

from sklearn.linear_model import MultiTaskElasticNetCV, MultiTaskElasticNet, ElasticNet, ElasticNetCV
from sklearn.feature_selection import SelectPercentile, f_regression
from sklearn.model_selection import cross_val_score
from sklearn import svm


import matplotlib.pyplot as plt


class FeatureSelection():
    def __init__(self):
        ##location info
        self.female_path = glob.glob("/home/hyunbin/seqff/output/female_notfiltered/normalizedbincount/*")
        self.male_path = glob.glob("/home/hyunbin/seqff/output/male_notfiltered/normalizedbincount/*")
        self.supple_path = "/home/hyunbin/seqff/supplementary-table2.csv"
        self.alluseablebins_path = "/home/hyunbin/seqff/output/male_notfiltered/alluseablebins"
        self.male_ff_path = "/home/hyunbin/seqff/output/ff_chrY.csv"

    def pearson_corr(self, all_auto_array):
        pearson_result = np.corrcoef(all_auto_array[:, 1000:])
        print("pearson result is completely calculated")

    def selection_with_corr(self):
        raise NotImplementedError

    def elastic_net_hyperparameter_opt(self, all_auto_array, male_ff_array):
        """
        tentative, so slow
        :param all_auto_array:
        :param male_ff_array:
        :return:
        """
        start = time.time()

        l1_ratio = 0.99  # could be changed
        regr = ElasticNetCV(l1_ratio=l1_ratio, n_alphas=100,
                                     cv=10)
        regr.fit(all_auto_array[:, 1000:].transpose(), male_ff_array[:, 1])
        end = time.time()

        elapsed = end - start

        print("alpha:", regr.alpha_)
        print("l1_ratio_:", regr.l1_ratio_)
        print(regr.coef_.shape)

        print("%d" % elapsed)

        return regr

    def selection_with_elastic_net(self, alpha, l1_ratio, all_auto_array, male_ff_array, save=0):
        regr = ElasticNet(alpha=alpha, l1_ratio=l1_ratio)
        regr.fit(all_auto_array[:,1000:].transpose(), male_ff_array[:, 1])

        print(regr.coef_)
        print(regr.coef_.shape)
        print(regr.intercept_)
        print(regr.intercept_.shape)

        #TEMP for saving file
        if save:
            with open("/home/hyunbin/seqff/output/feature_selection/enet_coef.csv", 'w') as outfile:
                for i in regr.coef_:
                    print(i, file=outfile)

    def selection_with_percentile(self, all_auto_array, male_ff_array, percentile=10):
        model = SelectPercentile(f_regression, percentile=percentile) #select 10% most significnat features
        model.fit(all_auto_array[:,1000:].transpose(), male_ff_array[:, 1])
        scores = -np.log10(model.pvalues_)
        scores /= scores.max()
        # plt.bar(X_indices - .45, scores, width=.2,
        #         label=r'Univariate score ($-Log(p_{value})$)', color='darkorange',
        #         edgecolor='black')

    def weights_of_svm(self, all_auto_array, male_ff_array, selector):
        X_indices = np.arange(all_auto_array[:, 1000:].transpose().shape[-1])

        clf = svm.SVC(kernel='linear')
        clf.fit(all_auto_array[:,1000:].transpose(), male_ff_array[:, 1])

        svm_weights = (clf.coef_ ** 2).sum(axis=0)
        svm_weights = svm_weights.max()

        clf_selected = svm.SVC(kernel='linear')
        clf_selected.fit(selector.transform())

        plt.bar(X_indices - .25, svm_weights, width=.2, label='SVM weight',
                color='navy', edgecolor='black')

        svm_weights_selected = (clf_selected.coef_ ** 2).sum(axis=0)
        svm_weights_selected = svm_weights_selected.max()

        plt.bar(X_indices[selector.get_support()] - .05, svm_weights_selected,
                width=.2, label='SVM weights after selection', color='c',
                edgecolor='black')

        plt.title("Comparing feature selection")
        plt.xlabel('Feature number')
        plt.yticks(())
        plt.axis('tight')
        plt.legend(loc='upper right')
        plt.show()

    def result(self, select):
        supple = pd.read_csv(self.supple_path)
        alluseablebins = pd.read_csv(self.alluseablebins_path, header=None)

        supple['alluseablebins'] = alluseablebins

        female_file_list = [".".join(filename.split("\\")[-1].split(".")[:-1]) for filename in self.female_path]
        male_file_list = [".".join(filename.split("\\")[-1].split(".")[:-1]) for filename in self.male_path]

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

        with open(self.male_ff_path, 'r')as infile:
            #male_ff_array = np.array(csv.reader(infile))

            male_ff_list = []
            for line in infile:
                line = line.strip()
                elem = line.split(",")
                elem[0] = ".".join(elem[0].split(".")[:-1])
                male_ff_list.append([elem[0], elem[1]])

            male_ff_array = np.array(male_ff_list)

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

        all_file_array = np.concatenate((np.array(female_file_list), np.array(male_file_list)), axis=None)
        all_auto_df.columns = all_file_array
        all_auto_df = all_auto_df.transpose()

        female_y_concat = pd.concat([df['normalizedbincount'] for df in female_y_df_list], axis=1)
        male_y_concat = pd.concat([df['normalizedbincount'] for df in male_y_df_list], axis=1)

        all_y_df = pd.concat([female_y_concat, male_y_concat], axis=1)
        all_y_df.index = male_y_df_list[0]['Unnamed: 0']

        all_file_array = np.concatenate((np.array(female_file_list), np.array(male_file_list)), axis=None)
        all_y_df.columns = all_file_array

        all_y_df = all_y_df.transpose()

        all_auto_array = np.array([all_auto_df[i].values for i in all_auto_df.columns])

        if select == 'pearson':
            corr = self.pearson_corr(all_auto_array)
            #TODO
            #save result to csvfile

        elif select == 'enet':
            #enet_cv = self.elastic_net_hyperparameter_opt(all_auto_array, male_ff_array)
            #alpha = float(enet_cv.alpha_)
            #l1_ratio = float(enet_cv.l1_ratio)
            enet = self.selection_with_elastic_net(alpha=0.007095309, l1_ratio=0.99,
                                                   all_auto_array=all_auto_array, male_ff_array=male_ff_array, save=True)


        elif select == '':
            raise NotImplementedError

        else:
            raise ValueError



if __name__ == "__main__":
    featureSelection = FeatureSelection()
    featureSelection.result(select='enet')

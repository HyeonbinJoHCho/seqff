import os
import numpy as np
import random
import time
import glob

import pandas as pd
from scipy.stats import ttest_ind


class ChrYTTest:
    def __init__(self, prefix):
        self.prefix = prefix
        self.out_dir_path = os.path.join(self.prefix, "t_test")
        self.male_data_dir_path = os.path.join(self.prefix, 'adult_male/newtemp')
        self.female_data_dir_path = os.path.join(self.prefix, 'female/newtemp')

        self.bininfo = pd.read_csv("/home/hyunbin/seqff/supplementary-table2.csv")
        self.bininfo.rename(columns={'Unnamed: 0': 'binName'}, inplace=True)
        self.bininfo = self.bininfo[self.bininfo['CHR'] == 'chrY']

        if not os.path.exists(self.out_dir_path):
            os.mkdir(self.out_dir_path)

    def load_bincounts(self, location, gender, suffix=''):
        bincounts_list = []
        files = glob.glob(location + "/*.newtemp")
        if gender == "female":
            files = random.sample(files, 10)
            file_names = [{'filename': os.path.basename(f)} for f in files]
            file_names = pd.DataFrame(file_names)
            file_names.to_csv(os.path.join(self.out_dir_path, 'female.info_%s.csv' % suffix), index=False)

        for file in files:
            bincounts = pd.read_csv(file, sep='\t', header=None)
            bincounts.columns = ['binName', 'counts']
            tmp_df = pd.merge(self.bininfo, bincounts, how='left', on='binName')
            tmp_df['counts'] = tmp_df[['counts']] .fillna(0)
            bincounts_list.append(tmp_df)
        return bincounts_list

    def t_test(self, suffix=''):
        male_bincounts_list = self.load_bincounts(self.male_data_dir_path, gender="male", suffix=suffix)
        female_bincounts_list = self.load_bincounts(self.female_data_dir_path, gender="female", suffix=suffix)

        male_bincounts_list = [df[df['CHR'] == 'chrY'] for df in male_bincounts_list]
        female_bincounts_list = [df[df['CHR'] == 'chrY'] for df in female_bincounts_list]

        result_list = []
        for bin_name in self.bininfo['binName'].values:
            x1 = [df[df['binName'] == bin_name]['counts'].values[0] for df in male_bincounts_list]
            x2 = [df[df['binName'] == bin_name]['counts'].values[0] for df in female_bincounts_list]
            result = ttest_ind(x1, x2, equal_var=False)
            result_list.append({
                'binName': bin_name,
                't-value': result.statistic,
                'p-value': result.pvalue
            })
        outfilename = os.path.join(self.out_dir_path, "t_test_%s.csv" % suffix)
        t_test_df = pd.DataFrame(result_list)
        t_test_df.to_csv(outfilename, index=False)
        return t_test_df

    def t_test_avg(self, num):
        t_test_df = None
        for i in range(num):
            tmp_t_test_df = self.t_test(suffix=str(i))
            if i == 0:
                t_test_df = tmp_t_test_df
            else:
                t_test_df['t-value'] = t_test_df['t-value'] + tmp_t_test_df['t-value']
                t_test_df['p-value'] = t_test_df['p-value'] + tmp_t_test_df['p-value']

        t_test_df['t-value'] = t_test_df['t-value'] / num
        t_test_df['p-value'] = t_test_df['p-value'] / num
        outfilename = os.path.join(self.out_dir_path, "t_test_avg_%s.csv" % num)
        t_test_df.to_csv(outfilename, index=False)


if __name__ == "__main__":
    prefix = "/home/hyunbin/seqff/output"
    ttest = ChrYTTest(prefix)
    # ttest.t_test()
    ttest.t_test_avg(100)

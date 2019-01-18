import time, os
import pandas as pd
import numpy as np
import math
from skmisc.loess import loess

import matplotlib.pyplot as plt
import temp_loessPredict

from loads import load_rdata, load_bininfo, load_sam, load_counts, load_bam
from creates import create_newtemp, create_normalizedbincount, create_alluseablebins, create_bincounts, create_results

#todo list
#1. code transplantation(?) from R line by line
#2. extract each variable in SupplementalFile1.RData   (current)
#3. code refinement
#4. multiprocessing or multi-threading?


class SeqFF:
    def __init__(self,
                 bininfo_loc="/home/hyunbin/seqff/supplementary-table2.csv",
                 input_loc="/home/woong/TEST/bam_mapped_male_set2/",
                 rdata="/home/hyunbin/seqff/supplementary-file1.rdata",
                 output_loc="/home/hyunbin/seqff/output/male_20190114/"):
        #file locations

        self.bininfodata_loc = bininfo_loc
        self.output_loc = output_loc
        self.input_loc = input_loc
        self.rdata = rdata

        self.newtemp_loc = os.path.join(self.output_loc, "newtemp")
        self.alluseablebins_loc = os.path.join(self.output_loc, "alluseablebins")
        self.normalizedbincount_loc = os.path.join(self.output_loc, "normalizedbincount")
        self.bincounts_loc = os.path.join(self.output_loc, "bincounts")
        self.results_loc = os.path.join(self.output_loc, "results")

        #rdata parameters load
        params_seqff = load_rdata(rdata)
        self.B = params_seqff['B']
        self.mu = params_seqff['mu']
        self.parameter_1 = params_seqff['parameter.1']
        self.parameter_2 = params_seqff['parameter.2']
        self.elnetbeta = params_seqff['elnetbeta']
        self.elnetintercept = params_seqff['elnetintercept']

        self.progress = len(os.listdir(self.input_loc))

    def prediction(self, bincounts, B, mu, parameter_1, parameter_2):
        """
        len(bincounts) :
        len(bincounts_2) :
        """
        bincounts[np.isnan(bincounts)] = 0

        pattern = 'chr[0-9]'
        bincounts_2 = bincounts[(bincounts.index.str.contains(pattern))]
        bincounts_3 = bincounts_2 - mu

        y_hat = np.dot(pd.concat([pd.Series([1]), bincounts_3]), B)
        sum_yhat = sum(y_hat[(pd.isna(y_hat)) == False])
        sum_bincounts = sum(bincounts_2[(pd.isna(bincounts_2)) == False])
        y_hat_rep = sum_yhat / sum_bincounts

        ff_hat = (y_hat_rep + parameter_1) * parameter_2

        return ff_hat

    def seqff(self):
        """
        Returns values of seqff, wrsc, enet score

        supplementary files can be downloaded below:
        https://obgyn.onlinelibrary.wiley.com/doi/abs/10.1002/pd.4615


        :param bininfoData: location of supplementary table2.csv file
        :type bininfoData: String

        :param inputData: directory path or file location of inputdata ( ".sam" or ".bam" or ".newtemp")
        :type inputData: String

        :param rdata: location of supplementary .rdata file
        :type rdata: String

        :param output_lod: where result files are(total 4 directories will be created)
        :type output_lod: String

        :return: None
        """

        start = time.time()

        # load bininfo
        bininfo = load_bininfo(self.bininfodata_loc)

        # load input files
        if os.path.isdir(self.input_loc):
            input_list = [self.input_loc + x for x in os.listdir(self.input_loc)]

        elif os.path.isfile(self.input_loc):
            input_list = [self.input_loc]

        else:
            raise FileNotFoundError("error occurred : inputData is not a Directory or File")

        for i, file in enumerate(input_list):
            filetype = file.split(".")[-1]
            # filetype : 'sam' or 'bam' or 'newtemp'
            if 'sam' in filetype:
                bincount = load_sam(file)

            elif 'newtemp' in filetype:
                bincount = load_counts(file)
                file = file.replace(".newtemp", "")  # TEMP .newtemp -> .bam

            elif 'bam' in filetype:
                bincount = load_bam(file)

            else:
                continue

            #CREATE newtemp file in "output_loc"/newtemp/
            create_newtemp(bincount, file, self.newtemp_loc)

            newtemp = pd.DataFrame.from_dict(bincount, orient='index')
            newtemp.reset_index(level=0, inplace=True)
            newtemp.rename(columns={'index': 'binName', 0: 'counts'}, inplace=True)

            temp_bininfo = bininfo.copy(deep=True)
            temp_bininfo = temp_bininfo.merge(newtemp, on='binName',
                                              how='left')  # missing value : NaN, not NA in pandas
            temp_bininfo['counts'] = temp_bininfo['counts'].fillna(0)

            temp_bininfo.sort_values(by='binorder', inplace=True)
            temp_bininfo.reset_index(drop=True)

            ####DATA PROCESSING #######################
            autosomebinsonly = []
            for index in range(61927):
                boolean = (temp_bininfo['FRS'][index] != 'NA') and \
                          (float(temp_bininfo['GC'][index]) > 0.316) and \
                          (temp_bininfo['CHR'][index] != 'chrX') and \
                          (temp_bininfo['CHR'][index] != 'chrY')
                autosomebinsonly.append(boolean)
            autosomebinsonly = pd.Series(autosomebinsonly)

            alluseablebins = []
            for index in range(61927):
                boolean = (temp_bininfo['FRS'][index] != "NA") and (float(temp_bininfo['GC'][index]) > 0.316)
                alluseablebins.append(boolean)
            alluseablebins = pd.Series(alluseablebins)

            #CREATE alluseablebins file in "output_loc"/alluseablebins
            #create_alluseablebins(alluseablebins, file, self.alluseablebins_loc)

            sum_counts = pd.Series(temp_bininfo['counts'])
            sum_counts = sum_counts[autosomebinsonly].sum(skipna=True)

            autoscaledtemp = pd.Series(temp_bininfo['counts'].loc[(autosomebinsonly)],
                                       copy=True) / sum_counts  # NA-related code removed
            allscaledtemp = pd.Series(temp_bininfo['counts'].loc[(alluseablebins)], copy=True) / sum_counts

            gc_index = {}
            cnt = 0
            for index, isauto in enumerate(autosomebinsonly):
                if isauto:
                    if temp_bininfo['GC'].iat[index] in gc_index:
                        gc_index[temp_bininfo['GC'].iat[index]].append(float(autoscaledtemp.iat[cnt]))
                        cnt += 1

                    else:
                        gc_index[temp_bininfo['GC'].iat[index]] = [float(autoscaledtemp.iat[cnt])]
                        cnt += 1

            key_list = []
            val_list = []
            for key, val in gc_index.items():
                key_list.append(key)
                val_list.append(np.median(val))

            loess_var = loess(key_list, val_list)  # default span : 0.75
            loess_var.fit()
            # y = loess.loess_prediction(newData, loessVar)
            # temp_loessPredict.loess_debugging(loessVar)

            ###prediction###
            loess_x = [float(gc) for index, gc in enumerate(temp_bininfo['GC']) if (alluseablebins[index])]
            # print(temp_bininfo['GC'])
            loess_fitted = loess_var.predict(loess_x)
            loess_fitted = list(loess_fitted.values)
            # print(loess_fitted)

            median_autoscaledtemp = np.median(autoscaledtemp)
            median_autoscaledtemp = float(median_autoscaledtemp)  # for fixed constant

            normalizedbincount = [(x + (median_autoscaledtemp - loess_fitted[index])) for index, x in
                                  enumerate(allscaledtemp)]

            #CREATE normalizedbincount in "output_loc"/normalizedbincount
            create_normalizedbincount(normalizedbincount, file, self.normalizedbincount_loc)

            bincounts = pd.Series(data=np.repeat(a=0.0, repeats=61927), index=temp_bininfo['binName'], dtype=np.float64)

            sum_normalizedbincount = sum([val for val in normalizedbincount if not math.isnan(val)])
            sum_normalizedbincount = float(sum_normalizedbincount)  # deep copy temporarily

            cnt = 0
            for index, x in enumerate(alluseablebins):
                if x == True:
                    data = (normalizedbincount[cnt] / sum_normalizedbincount) * len(normalizedbincount)
                    bincounts.iat[index] = data
                    cnt += 1

            #CREATE bincounts in "output_loc"/bincounts
            create_bincounts(bincounts, file, self.bincounts_loc)

            wrsc = self.prediction(bincounts, self.B, self.mu, self.parameter_1, self.parameter_2)
            enet = np.dot(bincounts, (self.elnetbeta)) + (self.elnetintercept)
            ff = (wrsc+enet) / 2

            result_lines = list()
            result_lines.append("SeqFF\tEnet\tWRSC")
            result_lines.append("{}\t{}\t{}".format(ff, enet, wrsc))

            #CREATE results of seqff (seqff paper result covered) in "output_loc"/results
            create_results(result_lines, file, self.results_loc)

            end = time.time()
            elapsed = end - start
            h = int(elapsed) // 3600
            m = (int(elapsed) - (h * 3600)) // 60
            s = (int(elapsed) % 60)
            print("elapsed time: %d hr %d min %d sec" % (h, m, s))
            print("elapsed :", elapsed)
            print("progress : {} / {}".format(i + 1, self.progress))


if __name__ == "__main__":
    seqff = SeqFF()
    seqff.seqff()
    print("completed")






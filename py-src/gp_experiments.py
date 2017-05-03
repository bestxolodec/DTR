from __future__ import print_function, division
import os
import sys
import pandas as pd
import scipy as sp
import numpy as np
from itertools import product
from timeit import default_timer as timer

import readline
import rpy2.robjects as ro
from itr import fit_and_predict
ro.r.source("/home/nbuser/DTR/src/clean_sources.R")


# get_ko_learning_value_and_predicted_dose = ro.globalenv['GetKOLearningValueAndPredictedDose']
# CHECK IF THIS GOES OK
# chen_1 = ro.globalenv['Scenario1Enriched']
# print(chen_1(100, 10, 0))


def generate_tuples(d, flables, slabels, sc):
    for i, j, k in product(*[range(o) for o in d.shape]):
        yield (sc, flables[i], slabels[j],  d[i, j, k])


class Experiment(object):

    def __init__(self, exp_params):
        self.n_cov = 10
        self.s_factors_percs = exp_params["s_factors_percs"]
        self.s_factors = sp.stats.norm.ppf(self.s_factors_percs)  # 50 factors evenly splitted
        self.n_repeats = exp_params["n_repeats"]
        self.granularity = exp_params["granularity"]
        self.n_train_list = exp_params["n_train_list"]
        self.n_test = exp_params["n_test"]
        assert len(self.n_train_list) > 0
        self.scenario = exp_params["scenario"]
        if "chen1" in scenario.lower():
            self.get_data = ro.globalenv['Scenario1Enriched']
            self.n_cov = 30
        elif "chen2" in scenario.lower():
            self.get_data = ro.globalenv['Scenario2Enriched']
        elif "chen4" in scenario.lower():
            self.get_data = ro.globalenv['Scenario4Enriched']
        else:
            raise "Uknown scenario: " + str(self.scenario)
        self.save_prefix = exp_params["save_prefix"]
        self.pred_value_func = ro.globalenv['PredValueGeneral']
        self.results = None

    def run(self):
        data = np.zeros((len(self.n_train_list), self.s_factors.size, self.n_repeats))
        for i, n_train in enumerate(self.n_train_list):
            print(n_train, end="\t")
            start = timer()
            for k in range(self.n_repeats):
                ko_train = get_data(n_train, self.n_cov, 777+k)
                ko_test = get_data(self.n_test, self.n_cov, 777+k)
                fit_params = {"verbose": False}
                # returns A, V, model, save only Values
                data[i, :, k] = fit_and_predict(ko_train, ko_test, self.granularity, self.s_factors,
                                                self.pred_value_func, fit_params)[1]
            print("elapsed {:.2f} min".format((timer() - start) / 60))
        self.results = data
        return self

    def write_to_file(self):
        save_fname = "{}_{}_rep".format(self.scenario, self.n_repeats)
        save_path = os.path.join(self.save_prefix, save_fname)
        print("Writing raw results to {}.......".format(save_path), end=" ")
        np.save(save_path, self.results)
        print("Success")
        df = pd.DataFrame.from_records(generate_tuples(self.results, self.n_train_list,
                                                       100 * self.s_factors_percs, self.scenario))
        df.columns = ["scenario", "sample_size", "s_factor", "value_f"]
        csv_save_path = save_path + ".csv"
        print("Writing csv results to {}.......".format(csv_save_path), end=" ")
        df.to_csv(csv_save_path, index=False)
        print("Success")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Experimets with GP')
    parser.add_argument('-s', '--scenario', action='store', default="chen1",
                        dest='scenario', help='Scenario to be run')
    parser.add_argument('--n_train_list', type=int, nargs='+', default=[50, 100, 200, 400, 800])
    parser.add_argument('--n_repeats', type=int, default=50)
    parser.add_argument('--granularity', type=int, default=50)
    parser.add_argument('--n_test', type=int, default=1000)
    parser.add_argument('--save_prefix', type=str, default="/home/nbuser/DTR/")
    parser.add_argument('--s_factors_percs', type=float, nargs='+', default=np.arange(.5, 1, .01))
    print(vars(parser.parse_args()))
    experiment = Experiment(vars(parser.parse_args()))
    experiment.run().write_to_file()

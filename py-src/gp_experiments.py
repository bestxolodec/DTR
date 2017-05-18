from __future__ import print_function, division
import os
import sys
import pandas as pd
import scipy as sp
import numpy as np
from itertools import product
from timeit import default_timer as timer
import logging
from pprint import pformat

import readline
import rpy2.robjects as ro
from itr import fit_and_predict
ro.r.source("/home/nbuser/DTR/src/clean_sources.R")

LOG_LEVEL = logging.WARNING
LOG_FILE = "/tmp/experiment.log"

logger = logging.getLogger()
logger.setLevel(LOG_LEVEL)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
fh = logging.FileHandler(LOG_FILE)
fh.setLevel(LOG_LEVEL)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setLevel(LOG_LEVEL)
ch.setFormatter(formatter)
logger.addHandler(ch)


def generate_tuples(d, flables, slabels, sc):
    for i, j, k in product(*[range(o) for o in d.shape]):
        yield (sc, flables[i], slabels[j],  d[i, j, k])


class Experiment(object):

    def __init__(self, exp_params):
        # self.n_cov = 10
        self.algo = exp_params["algorithm"]
        self.s_factors_percs = exp_params["s_factors_percs"]
        self.s_factors = sp.stats.norm.ppf(self.s_factors_percs)  # 50 factors evenly splitted
        self.n_repeats = exp_params["n_repeats"]
        self.granularity = exp_params["granularity"]
        self.n_train_list = exp_params["n_train_list"]
        self.n_test = exp_params["n_test"]
        assert len(self.n_train_list) > 0
        self.pred_value_func = ro.globalenv["PredValueGeneral"]
        self.fit_params = exp_params["fit_params"]
        self.scenario = exp_params["scenario"]
        self.results = None
        self.save_prefix = exp_params["save_prefix"]

    def _make_fun_gen_data_by_scenario(self):
        if "chen1" in self.scenario.lower(): return ro.globalenv["Scenario1Enriched"]
        elif "chen2" in self.scenario.lower(): return ro.globalenv["Scenario2Enriched"]
        elif "chen4" in self.scenario.lower(): return ro.globalenv["Scenario4Enriched"]
        elif "shvechikov1" in self.scenario.lower(): return ro.globalenv['GetDataForShvechikov.1']
        elif "shvechikov2" in self.scenario.lower(): return ro.globalenv['GetDataForShvechikov.2']
        else: raise "Unknown scenario: " + str(self.scenario)

    def _make_fun_fit_and_predict_by_algo(self):
        if "lcsl" in self.algo.lower():
            return lambda train, test: fit_and_predict(train, test, self.granularity, self.s_factors,
                                                       self.pred_value_func, self.fit_params)
        if "owl" in self.algo.lower():
            fun = ro.globalenv['GetKOLearningValueAndPredictedDose']
            return lambda train, test: [np.array(i) for i in fun(train, test)]

    def run(self):
        data = np.zeros((len(self.n_train_list), self.s_factors.size, self.n_repeats))
        get_data = self._make_fun_gen_data_by_scenario()
        fit_and_predict = self._make_fun_fit_and_predict_by_algo()
        for i, n_train in enumerate(self.n_train_list):
            start = timer()
            for k in range(self.n_repeats):
                train = get_data(n_train, 777 + k)   # n_of_samples, seed
                test = get_data(self.n_test, 777 + k)   # n_of_samples, seed
                # returns (A, V, ...); we save only Values
                data[i, :, k] = fit_and_predict(train, test)[1]
            logging.warning("{}\telapsed {:.2f} min".format(n_train, (timer() - start) / 60))
        self.results = data
        return self

    def make_fname(self):
        return "algo{}_reps{}_rests{}_norm{}_Xstand{}_Ystand{}".format(
            self.scenario, self.n_repeats, self.fit_params["n_restarts"], self.fit_params["normalize"],
            self.fit_params["standardize_X"], self.fit_params["standardize_Y"])


    def write_to_file(self):
        save_fname = "algo-{}_{}_reps-{}_rests-{}_norm-{}_Xstand-{}_Ystand-{}".format(
            self.algo, self.scenario, self.n_repeats, self.fit_params["n_restarts"],
            self.fit_params["normalize"], self.fit_params["standardize_X"], self.fit_params["standardize_Y"])
        save_path = os.path.join(self.save_prefix, save_fname)
        # np_save_path = save_path + ".npy"
        # logger.warning("Writing raw results to {}  .......".format(np_save_path))
        # np.save(np_save_path, self.results)
        # logger.warning("Success")
        data_tuples = generate_tuples(self.results, self.n_train_list, 100 * self.s_factors_percs, self.scenario)
        df = pd.DataFrame.from_records(data_tuples)
        df.columns = ["scenario", "sample_size", "s_factor", "value_f"]
        csv_save_path = save_path + ".csv"
        logger.warning("Writing csv results to {}  .......".format(csv_save_path))
        df.to_csv(csv_save_path, index=False)
        logger.warning("Success")
        return self


def parse_fit_params_arg(s):
    fit_params = {}
    for kv in s.strip().split(","):
        k, v = kv.strip().split(":")
        if "true" in v.lower():
            v = True
        elif "false" in v.lower():
            v = False
        else:
            v = int(v)
        fit_params[k] = v
    return fit_params


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Experimets with GP",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--scenario", action="store", default="chen1",
                        dest="scenario", help="Scenario to be run")
    parser.add_argument("--n_train_list", type=int, nargs="+", default=[50, 100, 200, 400, 800],
                        help="list of train sample sizes")
    parser.add_argument("--n_repeats", type=int, default=50, help="number of repeats of data sampling")
    parser.add_argument("--granularity", type=int, default=50, help="number of points in a treatment grid")
    parser.add_argument("--n_test", type=int, default=1000, help="number of samples in test data")
    parser.add_argument("--save_prefix", type=str, default="/home/nbuser/DTR/",
                        help="default path to save simulation results")
    parser.add_argument("--algorithm", type=str, default="LCSL",
                        help="algorithm to make predictions: one of OWL, LCSL")
    parser.add_argument("--s_factors_percs", type=float, nargs="+", default=np.arange(.5, 1, .01),
                        help="which percentiles to consider for variance penalty")
    parser.add_argument("--fit_params", type=str, default="",
                        help="String with 'key:value'")
    args_in_dict = vars(parser.parse_args())
    args_in_dict["fit_params"] = parse_fit_params_arg(args_in_dict["fit_params"])
    logger.warning(pformat(args_in_dict))
    experiment = Experiment(args_in_dict).run().write_to_file()

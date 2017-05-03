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


# CHECK IF THIS GOES OK
# chen_1 = ro.globalenv['Scenario1Enriched']
# print(chen_1(100, 10, 0))

def generate_tuples(d, flables, slabels, sc):
    for i, j, k in product(*[range(o) for o in d.shape]):
        yield (sc, flables[i], slabels[j],  d[i, j, k])


def get_params_for_scenario(scenario):
    n_cov, n_train_samples, n_test = 10, [50, 100, 200, 400, 800], 1000  # data gen
    s_factors_percs = np.arange(.5, 1, .01)
    s_factors = sp.stats.norm.ppf(s_factors_percs)  # 50 factors evenly splitted
    n_repeats = 50
    granularity = 50
    if "chen1" in scenario.lower():
        get_data = ro.globalenv['Scenario1Enriched']
        n_cov = 30
    elif "chen2" in scenario.lower():
        get_data = ro.globalenv['Scenario2Enriched']
    elif "chen4" in scenario.lower():
        get_data = ro.globalenv['Scenario4Enriched']
    else:
        raise "Uknown scenario"
    return get_data, n_cov, n_train_samples, n_test, s_factors_percs, s_factors, n_repeats, granularity


def run_simulation(scenario, save_prefix="/home/nbuser/DTR/"):
    params = get_params_for_scenario(scenario)
    get_data, n_cov, n_train_samples, n_test, s_factors_percs, s_factors, n_repeats, granularity = params
    data = np.zeros((len(n_train_samples), s_factors.size, n_repeats))
    for i, n_train in enumerate(n_train_samples):
        print(n_train, end="\t")
        start = timer()
        for k in range(n_repeats):
            ko_train = get_data(n_train, n_cov, 777+k)
            ko_test = get_data(n_test, n_cov, 777+k)
            data[i, :, k] = fit_and_predict(ko_train, ko_test, granularity,
                                            s_factors, pred_value_general)[1]
        print("elapsed {:.2f} min".format((timer() - start) / 60))
    save_fname = "{}_{}_rep".format(scenario, n_repeats)
    save_path = os.path.join(save_prefix, save_fname)
    np.save(save_path, data)
    df_data = pd.DataFrame.from_records(generate_tuples(data, n_train_samples, s_factors_percs, scenario))
    df_data.to_csv(save_path + ".csv", index=False)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Experimets with GP')
    parser.add_argument('-s', '--scenario', action='store',
                        dest='scenario', help='Scenario to be run')
    args = parser.parse_args()
    run_simulation(args.scenario)

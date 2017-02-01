import GPy as gpy
from GPy.mappings import Linear, Constant, Additive
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from collections import Iterable

def get_initial_inducing(n_of_inducing, X):
    from sklearn.cluster import KMeans
    return KMeans(n_clusters=n_of_inducing, random_state=0, copy_x=False,
                  max_iter=100, n_init=2).fit(X).cluster_centers_

    assert X.ndim == 2
    results = []
    for i in range(X.shape[1]):
        results.append(np.random.uniform(*np.percentile(X[:, i], [0, 100]), size=n_of_inducing))
    return np.vstack(results).T


def unpack_from_rpyobj(obj):
    return [np.array(obj.rx2(name)) for name in ['covariates', 'treatment', 'reward', "optimal.treatment"]]


def unpack_data_from_rpyobj(obj):
    if 'covariates' in obj.names:
        unpacked_seq = [np.array(obj.rx2(name)) for name in ['covariates', 'treatment', 'reward', "optimal.treatment"]]
    else:
        unpacked_seq = [np.array(obj.rx2(name)) for name in ['X', 'A', 'R', 'D_opt']]
    return [ar.reshape(-1, 1) if ar.ndim == 1 else ar for ar in unpacked_seq]


def fit_GP(X, Y, fit_params):
    # :fit_params: dict with – mean_fn, n_restarts, n_inducing, verbose, robust
    kern = gpy.kern.RBF(X.shape[1],  ARD=True)  # n_of_dimensions
    mf = Additive(Linear(X.shape[1],Y.shape[1]), Constant(X.shape[1],Y.shape[1]))
    n_inducing = fit_params.get("n_inducing")
    if n_inducing is not None:  # doing sparse regression
        Z = get_initial_inducing(n_inducing, X)
        # TODO: mean_function is not supported in default constructor
        m = gpy.models.SparseGPRegression(X, Y, kernel=kern, Z=Z)   # , mean_function=mf)
    else:
        mf = mf if fit_params.get("mean_fn") is not None else None
        m = gpy.models.GPRegression(X, Y, kern, mean_function=mf)
    # m.optimize()  # TODO: optimization  with restarts
    n_restarts, robust, verbose = [fit_params.get(i) for i in ["robust", "verbose"]]
    n_restarts = fit_params.get("n_restarts", 1)
    m.optimize_restarts(num_restarts=n_restarts, num_processes=min(4, n_restarts), robust=robust, verbose=verbose)
    return m

def fit_GP(X, Y, n_restarts=1, n_inducing=None):
    kern = gpy.kern.RBF(X.shape[1],  ARD=True)  # n_of_dimensions
    mf = Additive(Linear(X.shape[1],Y.shape[1]), Constant(X.shape[1],Y.shape[1]))
    if n_inducing is not None:  # doing sparse regression
        Z = get_initial_inducing(n_inducing, X)
        # TODO: mean_function is not supported in default constructor
        m = gpy.models.SparseGPRegression(X, Y, kernel=kern, Z=Z)   # , mean_function=mf)
    else:
        m = gpy.models.GPRegression(X, Y, kern, mean_function=mf)
    # m.optimize()  # TODO: optimization  with restarts
    m.optimize_restarts(num_restarts=n_restarts, num_processes=min(4, n_restarts), robust=True, verbose=True)
    return m


def predict_with_GP_lower_surface(m, X_test, s_vec, granularity, space):
    assert isinstance(s_vec, Iterable)  # s_vec should be at least python list or np.array
    prediction = []
    n_in_batch = 100 * granularity  # per batch prediction for memory saving
    split_indexes = np.arange(n_in_batch, X_test.shape[0], n_in_batch)
    for batch in np.array_split(X_test, split_indexes, axis=0):
        ms, vs = [o.reshape(-1, granularity) for o in m.predict(batch)]
        # lower_surf.shape = (s_factors, num of test objects, num of grid points)
        lower_surf = ms - np.asarray(s_vec).reshape(-1, 1, 1) * np.sqrt(vs)
        prediction.append(space[lower_surf.argmax(axis=-1).T])
    return np.vstack(prediction) # shape = (num of test objects, s_factors)


def get_mesh_of_cov_treat(C, a_min, a_max, granularity):
    space = np.linspace(a_min,a_max,granularity)
    n_test = C.shape[0]
    return np.hstack([np.repeat(C, granularity, axis=0), np.tile(space, n_test).reshape(-1,1)]), space


def get_brute_treatment_prediction(granularity, a_min, a_max, X, R, C_test, s_vec, n_restarts=1, n_inducing=None):
    X_mesh, space = get_mesh_of_cov_treat(C_test, a_min, a_max, granularity)
    m = fit_GP(X, R, n_restarts=n_restarts,  n_inducing=n_inducing)
    return predict_with_GP_lower_surface(m, X_mesh, s_vec, granularity, space), m


def np2r(array):
    from rpy2.robjects import r
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    nr, nc = array.shape if array.ndim > 1 else (array.size, 1)
    array_r = r.matrix(array, nrow=nr, ncol=nc)
    return r.assign("array_r", array_r)


def fit_and_predict(train, test, granularity, s_vec, get_prediction, n_restarts=1, n_inducing=None, eps=0.05):
    # :param: get prediction  – function (prediction, test) which returns value for prediction
    C, A, R, A_opt = unpack_data_from_rpyobj(train)
    C_test, A_test, R_test, A_opt_test = unpack_data_from_rpyobj(test)
    X, X_test = [np.hstack([a, b]) for a, b in zip([C, C_test], [A, A_test])]
    a_min, a_max = np.percentile(A, [0, 100]) + [-eps, +eps]
    A_preds, m = get_brute_treatment_prediction(
        granularity, a_min, a_max, X, R, C_test, s_vec, n_restarts=n_restarts, n_inducing=n_inducing)
    # for each column (treatment predicted for particular s) get value
    values = [get_prediction(np2r(A_pred.reshape(-1,1)), test)[0] for A_pred in A_preds.T]
    return A_preds, np.array(values), m


def get_kc_exp_covar(ls_c, c_t, C, sigma_2f):
    # length-scale, new covariate, train covariates, variance of f
    return sigma_2f * np.exp(- .5 * np.square(c_t - C) / ls_c).reshape(-1)


def rbf_outer(x_row, x_col, ls, beta=1):
    return np.exp(- .5 * np.square(x_row.reshape(-1, 1) - np.array(x_col).ravel()) / ls * beta)


def mean_profile(space, A, ls_a, alphas, beta):
    return np.sum(rbf_outer(A, space, ls_a, beta) * alphas, axis=0)


def var_profile(space, A, ls_a, gammas, sigma_2f):
    o = np.square(space - A).T
    idx = range(len(space))
    exp_space_A_A = np.exp(- .5 * np.add.outer(o, o)[idx, :, idx, :] / ls_a)
    return sigma_2f - (exp_space_A_A * gammas).sum(axis=(-2, -1))


def grad_mean_profile(x, A, ls_a, alphas, beta):
    x = np.array(x)
    assert x.size == 1
    return (alphas * (A.reshape(-1, 1) - x) / ls_a * beta * rbf_outer(A, x, ls_a, beta)).sum()


def plot_curve_on_space(space, A, ls_a, alphas, beta, a_max, a_min=0):
    mean_profile_values = mean_profile(space, A, ls_a, alphas, beta)
    plt.plot(space, mean_profile_values)
    plt.title("Argmax = {}".format(space[np.argmax(mean_profile_values)]))


def generate_sample(n_of_train, a_max=100, a_min=0, verbose=False, seed=0):
    if seed is not None:
        np.random.seed(seed=seed)
    space = np.linspace(a_min, a_max, 10000)
    A = np.random.uniform(a_min, a_max, size=n_of_train)
    alphas = sp.stats.norm(0, 7).rvs(n_of_train).reshape(-1, 1)
    ls_a = sp.stats.gamma(1, scale=3).rvs(1)
    return space, A, alphas, ls_a


print("itr.py is imported!")

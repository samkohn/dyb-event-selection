'''
Rate-only analysis.

'''
import math

def n_obs_AD(ibd_rate, livetime):
    return ibd_rate * livetime

def n_obs_EH(ibd_rate_dict, livetime_dict):
    n_tot = 0
    for det in ibd_rate_dict:
        n_tot += n_obs_AD(ibd_rate_dict[det], livetime_dict[det])
    return n_tot

def n_obs_EH_error(ibdrate_errors_dict, livetime_dict):
    error = 0
    for det in ibdrate_errors_dict:
        error += pow(ibdrate_errors_dict[det] * livetime_dict[det], 2)
    error = math.sqrt(error)
    return error

def n_far_exp(w, n_EH1_obs, n_EH2_obs):
    return w['EH1']*n_EH1_obs + w['EH2']*n_EH2_obs

def n_far_exp_error(w, EH1_error, EH2_error):
    return math.sqrt(n_far_exp(w, EH1_error, EH2_error))

def r_measured(ibd_rates, livetimes, ibd_rate_errors):
    '''
    ibd_rates = {'EH1': {1: r1, 2: r2}, 'EH2': ...}

    livetimes = {'EH1': {1: t1, 2: t2}, 'EH2': ...}

    '''
    n_obs = {EH: n_obs_EH(ibd_rates[EH], livetimes[EH]) for EH in ibd_rates}
    n_obs_errors = {EH: n_obs_EH_error(ibd_rate_errors[EH], livetimes[EH]) for
            EH in ibd_rate_errors}
    fkl = create_fkl_table(livetimes)
    w = create_w_weights(fkl)
    exp_EH3 = n_far_exp(w, n_obs['EH1'], n_obs['EH2'])
    exp_EH3_error = n_far_exp_error(w, n_obs_errors['EH1'],
            n_obs_errors['EH2'])
    obs_EH3 = n_obs['EH3']
    obs_EH3_error = n_obs_errors['EH3']
    r = obs_EH3/exp_EH3
    error = math.sqrt(pow(obs_EH3_error/obs_EH3, 2) +
            pow(exp_EH3_error/exp_EH3, 2))
    return r, error


distances = {
        'D': [{
            'EH1': [362.38, 357.94],
            'EH2': [1332.48, 1337.43],
            'EH3': [1919.63, 1917.52, 1925.26, 1923.15]
            }, {
            'EH1': [371.76, 368.41],
            'EH2': [1358.15, 1362.88],
            'EH3': [1894.34, 1891.98, 1899.86, 1897.51]
            }],
        'L': [{
            'EH1': [903.47, 903.35],
            'EH2': [467.57, 472.97],
            'EH3': [1533.18, 1534.92, 1538.93, 1540.67]
            }, {
            'EH1': [817.16, 816.90],
            'EH2': [489.58, 495.35],
            'EH3': [1533.63, 1535.03, 1539.47, 1540.87]
            }, {
            'EH1': [1353.62, 1354.23],
            'EH2': [557.58, 558.71],
            'EH3': [1551.38, 1554.77, 1556.34, 1559.72]
            }, {
            'EH1': [1265.32, 1265.89],
            'EH2': [499.21, 501.07],
            'EH3': [1524.94, 1528.05, 1530.08, 1533.18]
            }]
        }

def create_fkl_table(livetimes, from_paper=False):
    if from_paper:
        return {
            'EH1': {'D':3.5022, 'L':0.9255},
            'EH2': {'D':0.2338, 'L':3.4333},
            'EH3': {'D':0.2423, 'L':0.7577}
            }
    fkl = {
            'EH1': {'D': 0, 'L': 0},
            'EH2': {'D': 0, 'L': 0},
            'EH3': {'D': 0, 'L': 0}
            }
    for EH in fkl:
        for key in distances:
            for distance_map in distances[key]:
                listing = distance_map[EH]
                for i, distance in enumerate(listing):
                    fkl[EH][key] += pow(distance, -2) * livetimes[EH][i+1]

    EH3_tot = fkl['EH3']['D'] + fkl['EH3']['L']
    for EH in fkl:
        for reactor in fkl[EH]:
            fkl[EH][reactor] /= EH3_tot
    return fkl

def create_w_weights(fkl_table):
    fkl = fkl_table
    denominator = (fkl['EH1']['D']*fkl['EH2']['L'] -
            fkl['EH1']['L']*fkl['EH2']['D'])
    w_1 = ((fkl['EH3']['D']*fkl['EH2']['L'] - fkl['EH3']['L']*fkl['EH2']['D'])
            /denominator)
    w_2 = ((fkl['EH3']['L']*fkl['EH1']['D'] - fkl['EH3']['D']*fkl['EH1']['L'])
            /denominator)
    return {'EH1': w_1, 'EH2': w_2}

def create_beta_weights(fkl_table):
    fkl = fkl_table
    beta_denom = fkl['EH3']['D'] + fkl['EH3']['L']
    beta_1 = (fkl['EH1']['D'] + fkl['EH1']['L'])/beta_denom
    beta_2 = (fkl['EH2']['D'] + fkl['EH2']['L'])/beta_denom
    return {'EH1': beta_1, 'EH2': beta_2}

eta = {'EH1': 0.18, 'EH2': 0.206, 'EH3': 0.789}

def sin22theta13(ibd_rates, livetimes, ibd_rate_errors):
    r, r_error = r_measured(ibd_rates, livetimes, ibd_rate_errors)
    fkl = create_fkl_table(livetimes)
    w = create_w_weights(fkl)
    beta = create_beta_weights(fkl)
    eta_far = eta['EH3']
    eta_near = w['EH1']*beta['EH1']*eta['EH1'] +w['EH2']*beta['EH2']*eta['EH2']
    theta13 = (1-r)/(eta_far - r*eta_near)
    theta13_error = math.sqrt(2) * r_error
    return (theta13, theta13_error)


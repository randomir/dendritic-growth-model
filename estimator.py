#!/usr/bin/env python
# -*- python-version >= 3.4 -*-
from collections import OrderedDict

import numpy as np
from scipy.optimize import minimize

from model import InvalidModelParam
from stats import neuronset_paths, neuronset_basal_stats
from main import simulate


if __name__ == '__main__':
    young_stats = neuronset_basal_stats(neuronset_paths['young'])
    old_stats = neuronset_basal_stats(neuronset_paths['old'])

    all_params = OrderedDict([
        ('B', 1), ('E', 1), ('S', 0), ('N_be', 100), ('N_e', 0),
        ('offset_in', 0), ('mean_in', 1), ('sd_in', 1),
        ('offset_be', 0), ('mean_be', 1), ('sd_be', 1),
        ('offset_e', 0), ('mean_e', 1), ('sd_e', 1)])
    
    # all relevant params at once
    # optim_params = (
    #     'B', 'E', 'S', 
    #     'offset_in', 'mean_in', 'sd_in',
    #     'offset_be', 'mean_be', 'sd_be')
    # p0 = np.array([1, 1, 0,
    #                0, 1, 1,
    #                0, 1, 1])

    # only structural params
    optim_params = (
        'B', 'E', 'S', 
    )
    p0 = np.array([10, 0.5, 0])
    bounds = [(0, 100), (0, 1), (-1, 1)]

    def stats_error(a, b):
        return (
            (a['degree']['mean'] - b['degree']['mean'])**2 +
            (a['degree']['stdev'] - b['degree']['stdev'])**2 +
            (a['mean_order']['mean'] - b['mean_order']['mean'])**2 +
            (a['mean_order']['stdev'] - b['mean_order']['stdev'])**2 +
            (a['asymmetry_index']['mean'] - b['asymmetry_index']['mean'])**2 +
            (a['asymmetry_index']['stdev'] - b['asymmetry_index']['stdev'])**2 #+
            # (a['total_length']['mean'] - b['total_length']['mean'])**2 +
            # (a['total_length']['stdev'] - b['total_length']['stdev'])**2 +
            # (a['intermediate_lengths']['mean'] - b['intermediate_lengths']['mean'])**2 +
            # (a['intermediate_lengths']['stdev'] - b['intermediate_lengths']['stdev'])**2 +
            # (a['terminal_lengths']['mean'] - b['terminal_lengths']['mean'])**2 +
            # (a['terminal_lengths']['stdev'] - b['terminal_lengths']['stdev'])**2
        )
            
    def fn(p):
        all_params.update(zip(optim_params, p))
        try:
            p_stats = simulate(all_params, n=1000)
            return stats_error(p_stats, old_stats)
        except InvalidModelParam:
            return np.inf
    
    res = minimize(fn, p0, method='nelder-mead', options={'xtol': 1e-8, 'maxiter': 1000, 'disp': True})
    print(res.x)

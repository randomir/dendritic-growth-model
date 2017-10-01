#!/usr/bin/env python
# -*- python-version >= 3.4 -*-
from pprint import pprint

from utils import counted
from model import Segment, DendriticTree
from stats import map_with_stats


def simulate_and_measure(params):
    tree = DendriticTree(**params)
    tree.grow()
    return tree.stats()


def simulate(params, n):
    return map_with_stats(simulate_and_measure, [[params]] * n)


# testing

def run_single(params):
    tree = DendriticTree(**params)
    tree.grow()

    print(tree.root.pformat())
    print("Degree at root:", tree.root.degree)
    print("Tree asymmetry index:", tree.asymmetry_index)
    print("Total length:", tree.total_length)

    #print("Function calls", counted.called)
    #print("Function times", counted.timing)


def run_multi(params, n):
    stats = simulate(params, n)
    print("Tree simulation stats (%d trees):" % n)
    pprint(stats)


if __name__ == '__main__':
    # S1-Rat Cortical Layer 2/3 Pyramidal Cell Basal Dendrites
    pyramidal_params = dict(
        B=2.52, E=0.73, S=0.5, N_be=312, N_e=96,
        offset_in=0, mean_in=6, sd_in=5,
        offset_be=0, mean_be=0.2, sd_be=0.2*0.47,
        offset_e=0, mean_e=0.86, sd_e=0.86*0.47
    )

    # Guinea Pig Purkinje Cell Dendritic Tree
    purkinje_params = dict(
        B=95, E=0.69, S=-0.14, N_be=10, N_e=0,
        offset_in=0.7, mean_in=10.63, sd_in=7.53
    )

    run_single(pyramidal_params)
    run_multi(pyramidal_params, 1000)

    # run_single(purkinje_params)
    # run_multi(purkinje_params, 1000)

#!/usr/bin/env python
# -*- python-version >= 3.4 -*-
"""
For a set of SWC neurons (subset of dendrites), load them in the model [1] and 
calculate stats which can later be used for parameter estimation.
"""
import statistics
import glob
from pprint import pprint
from plucky import merge

from swc import LineSegment, read_neuron
from model import Segment, DendriticTree


def map_with_stats(fn, argset, verbose=False):
    """Run `fn` over a list of parameters (positional arguments) in `argset`,
    calculating stats of all values present in response of each run.
    """
    def add_listified(x, y):
        if not isinstance(x, list):
            x = [x]
        if not isinstance(y, list):
            y = [y]
        return x + y

    sums = {}
    for args in argset:
        measures = fn(*args)
        if verbose:
            print("Run with args %r measures:" % args, measures)
        sums = merge(sums, measures, add_listified, recurse_list=False)

    stats = {}
    for k, v in sums.items():
        stats[k] = dict(total=sum(v),
                        mean=statistics.mean(v),
                        median=statistics.median(v),
                        stdev=statistics.stdev(v))
    return stats


def get_apical_linesegments(neuron):
    # get the one and only apical dendrite
    for c in neuron.root.children:
        if c.type == LineSegment.APICAL:
            return c

def get_basal_linesegments_set(neuron):
    # get all basal dendrites generator (of root line segments)
    for c in neuron.root.children:
        if c.type == LineSegment.BASAL:
            yield c


def build_dendrite_from_linesegments(root_linesegment):
    dendrite = DendriticTree()
    dendrite.empty()
    
    def trace(ls, parent_segment, parent_order, length=0):
        length += ls.length
        n_children = len(ls.children)
        if n_children == 0:
            # this is a terminal segment
            segment = Segment(dendrite, parent_order+1, parent_segment)
            # TODO: how to split total length we have into initial and elongated?
            segment.initial_len = length
            dendrite.terminal_segments.add(segment)
            return segment
        elif n_children == 1:
            # intermediate line-segment without branches (invalid in our model),
            # is still a segment growing...
            return trace(ls.children[0], parent_segment, parent_order, length)
        elif n_children == 2:
            # branching; finish tracing this segment and fork
            # (this is an intermediate segment)
            segment = Segment(dendrite, parent_order+1, parent_segment)
            segment.initial_len = length
            segment.children = [
                trace(ls.children[0], segment, segment.order, 0),
                trace(ls.children[1], segment, segment.order, 0)
            ]
            dendrite.intermediate_segments.add(segment)
            return segment
        else:
            raise Exception("Invalid LineSegment tree (3-way branch)")
    
    dendrite.root = trace(root_linesegment, parent_segment=None, parent_order=0)

    for segment in dendrite.terminal_segments:
        segment.update_degree()

    return dendrite



def load_dendrite_from_swc(path):
    neuron = read_neuron(path)
    apical = get_apical_linesegments(neuron)
    dendrite = build_dendrite_from_linesegments(apical)
    return dendrite


def test_load_dendrite():
    dendrite = load_dendrite_from_swc('../data/smit-rigter-mouse/92-1631.CNG.swc')
    print(dendrite.root.pformat())
    print("Degree at root:", dendrite.root.degree)
    print("Tree asymmetry index:", dendrite.asymmetry_index)
    print("Total length:", dendrite.total_length)
    print("Stats:", dendrite.stats())



def apical_dendrite_stats_from_swc(path):
    neuron = read_neuron(path)
    apical = get_apical_linesegments(neuron)
    dendrite = build_dendrite_from_linesegments(apical)
    return dendrite.stats()

def youngest_neurons_apical_stats():
    # youngest neurons 92-* (9 days, 15 neurons)
    paths = glob.glob('../data/smit-rigter-mouse/ws*.CNG.swc')
    stats = map_with_stats(apical_dendrite_stats_from_swc, [[path] for path in paths])
    print("Apical dendrites (9-day old) stats (%d neurons):" % len(paths))
    pprint(stats)

def oldest_neurons_apical_stats():
    # oldest neurons 92-* (365 days, 19 neurons)
    paths = glob.glob('../data/smit-rigter-mouse/92-*.CNG.swc')
    stats = map_with_stats(apical_dendrite_stats_from_swc, [[path] for path in paths])
    print("Apical dendrites (365-day old) stats (%d neurons):" % len(paths))
    pprint(stats)



def basal_dendrites_iter(paths):
    for path in paths:
        neuron = read_neuron(path)
        basal_linesegments = get_basal_linesegments_set(neuron)
        for ls in basal_linesegments:
            dendrite = build_dendrite_from_linesegments(ls)
            yield dendrite

def neuronset_basal_stats(paths):
    dendrites_argset = [[d] for d in basal_dendrites_iter(paths)]
    stats = map_with_stats(lambda d: d.stats(), dendrites_argset)
    print("Basal dendrites stats ({} neurons, {} dendrites):".format(len(paths), len(dendrites_argset)))
    pprint(stats)


if __name__ == '__main__':
    # youngest neurons 92-* (9 days, 15 neurons)
    youngest_paths = glob.glob('../data/smit-rigter-mouse/ws*.CNG.swc')

    # oldest neurons 92-* (365 days, 19 neurons)
    oldest_paths = glob.glob('../data/smit-rigter-mouse/92-*.CNG.swc')

    youngest_neurons_apical_stats()
    oldest_neurons_apical_stats()

    neuronset_basal_stats(youngest_paths)
    neuronset_basal_stats(oldest_paths)
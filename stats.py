#!/usr/bin/env python
# -*- python-version >= 3.4 -*-
"""
For a set of SWC neurons (subset of dendrites), load them in the model [1] and 
calculate stats which can later be used for parameter estimation.
"""
import statistics
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


def get_apical(neuron):
    # get the first apical dendrite
    for c in neuron.root.children:
        if c.type == LineSegment.APICAL:
            return c


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


def load_dendrite_from_swc(filename):
    neuron = read_neuron(filename)
    apical = get_apical(neuron)
    dendrite = build_dendrite_from_linesegments(apical)
    return dendrite


def test_load_dendrite():
    dendrite = load_dendrite_from_swc('../data/smit-rigter-mouse/92-1631.CNG.swc')
    print(dendrite.root.pformat())
    print("Degree at root:", dendrite.root.degree)
    print("Tree asymmetry index:", dendrite.asymmetry_index)
    print("Total length:", dendrite.total_length)

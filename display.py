#!/usr/bin/env python

"""
Parse and plot neurons from [SWC file
format](http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html).

In short, SWC is a text format, with one segment per line, with comment line
starting with `#`.

Each segment-describing line has 7 space separated fields:
 - id
 - structure id (0: undefined, 1: soma, 2: axon, 3: (basal) dendrite, 4: apical dendrite, 5+: custom); usually used 1, 3 and 4
 - x [um]
 - y [um]
 - z [um]
 - radius [um]
 - parent id (-1 for root; parent always appears in file before any of its children)
"""

import re
import os
import sys
import logging
from collections import namedtuple

import matplotlib.pyplot as plt
import matplotlib.lines as mlines


logger = logging.getLogger('display')

SegmentGeometry = namedtuple('SegmentGeometry', 'x y z r')


class Segment(object):
    "Simple segment definition."

    # segment types (ids and names)
    SOMA = 1
    AXON = 2
    BASAL = 3
    APICAL = 4
    typename = {SOMA: 'Soma', AXON: 'Axon',
                BASAL: 'Basal dendrite', APICAL: 'Apical dendrite'}

    def __init__(self, id, type, geometry, parent=None, children=None):
        self.id = id
        self.type = type
        self.geom = geometry
        self.parent = parent
        self.children = children or set()


class Neuron(object):
    name = None
    root = None
    # ref map: segment id -> segment object, for easier lookup/construction
    segments = dict()

    def add_segment_from_raw(self, id, type, x, y, z, r, parent_id):
        geom = SegmentGeometry(x, y, z, r)
        if self.root is None:
            if parent_id == -1 and type == 1:
                self.segments[id] = self.root = Segment(id, type, geom)
                return self.root
            else:
                raise Exception("Expecting soma root segment first.")
        
        try:
            parent = self.segments[parent_id]
        except:
            raise Exception("Invalid reference to parent segment.")
        
        self.segments[id] = segment = Segment(id, type, geom, parent)
        parent.children.add(segment)
        return segment


def read_neuron(filename):
    neuron = Neuron()
    neuron.name = os.path.basename(filename)
    with open(filename, 'r') as fin:
        for line in fin:
            if not re.match('^\d+', line):
                continue
            fields = line.split()
            if len(fields) != 7:
                logger.warning("Invalid line %r while parsing %r.", line, filename)
                continue
            fns = (int, int, float, float, float, float, int)
            values = [f(v) for f, v in zip(fns, fields)]
            neuron.add_segment_from_raw(*values)
    return neuron


def draw_neuron(neuron):
    segment_color = {Segment.SOMA: 'r', Segment.AXON: 'k',
                     Segment.BASAL: 'b', Segment.APICAL: 'g'}

    fig, ax = plt.subplots()
    # total dendritic length
    tdl = 0

    # plot each segment with line from parent to segment
    for segment in neuron.segments.values():
        parent = segment.parent
        if not parent:
            continue
        start = parent.geom
        end = segment.geom
        ax.annotate(str(segment.id), xy=(end.x, end.y), color='#cccccc')
        ax.plot([start.x, end.x], [start.y, end.y], 'o--',
                color=segment_color[segment.type])
        if segment.type in (Segment.BASAL, Segment.APICAL):
            tdl += ((end.x - start.x)**2 + (end.y - start.y)**2 + (end.z - start.z)**2)**0.5

    # create legend lines and labels
    legend_handles = []
    for typ in sorted(segment_color.keys()):
        legend_handles.append(
            mlines.Line2D(
                [], [], marker='o', linestyle='--',
                color=segment_color[typ], label=Segment.typename[typ]))
    plt.legend(handles=legend_handles)

    plt.title("%s (total dendritic length = %d)" % (neuron.name, tdl))

    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: %s NEURON.SWC" % __file__)
        sys.exit(1)
    filename = sys.argv[1]
    neuron = read_neuron(filename)
    draw_neuron(neuron)

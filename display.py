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
import logging
from collections import namedtuple

import matplotlib.pyplot as plt


logger = logging.getLogger('display')

SegmentGeometry = namedtuple('SegmentGeometry', 'x y z r')


class Segment(object):
    "Simple segment definition."

    def __init__(self, id, type, geometry, parent=None, children=None):
        self.id = id
        self.type = type
        self.geom = geometry
        self.parent = parent
        self.children = children or set()


class Neuron(object):
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


#neuron = read_neuron('../data/smit-rigter-mouse/92-1521.CNG.swc')

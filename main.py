"""
Dendritic Geometry Model, described by Jaap van Pelt et al. in [1], chapter 7.1.

[1] Computational neuroscience: Realistic modeling for experimentalists,
edited by Erik De Schutter, 2001, CRC Press.
"""
# -*- python-version >= 3.3 -*-
import math
import random
import textwrap

from utils import counted


class Segment(object):
    # segment's level in the dendritic tree (gamma param in [1])
    centrifugal_order = 0
    
    # segment's degree: the number of terminal segments in a subtree rooted at this segment.
    # if segment is a terminal segment, it has a degree of one
    degree = 1

    # link to root / dendrite tree
    # (holds growth params and current tree properties, like number of terminal segments)
    dendrite = None
    
    # link to parent, for easier traversal
    parent = None
    
    # list of daughter segments
    children = None

    initial_len = 0
    elongated_len = 0
    @property
    def total_len(self):
        return self.initial_len + self.elongated_len
    
    @counted
    def __init__(self, dendrite, order=0, parent=None):
        self.dendrite = dendrite
        self.parent = parent
        self.centrifugal_order = order

    def branching_probability_normalization_constant(self):
        s = 0
        for terminal in self.dendrite.terminal_segments:
            s += math.pow(2, -self.dendrite.parameters['S'] * terminal.centrifugal_order)
        return len(self.dendrite.terminal_segments) / s

    @counted
    def branch(self):
        if self.children:
            return
        
        C = self.branching_probability_normalization_constant()
        B = self.dendrite.parameters['B']
        E = self.dendrite.parameters['E']
        S = self.dendrite.parameters['S']
        N = self.dendrite.parameters['N']
        n_i = len(self.dendrite.terminal_segments)
        gamma = self.centrifugal_order
        p_i = C * math.pow(2, -S * gamma) * B / N / math.pow(n_i, E)
        
        if random.random() > p_i:
            # no branching
            return
        
        self.children = [Segment(self.dendrite, self.centrifugal_order + 1, self),
                         Segment(self.dendrite, self.centrifugal_order + 1, self)]
        self.update_degree()
        self.dendrite.terminal_segments.remove(self)
        self.dendrite.terminal_segments.update(self.children)

    @counted
    def update_degree(self):
        if self.children:
            self.degree = sum([c.degree for c in self.children])
        else:
            self.degree = 1
        if self.parent:
            self.parent.update_degree()

    def pformat(self):
        if self.children:
            children = textwrap.indent(",\n".join([c.pformat() for c in self.children]), prefix=" "*4)
            children = "\n%s\n" % children
        else:
            children = ""
        return "Segment(order=%s, degree=%s, children=[%s])" % (self.centrifugal_order, self.degree, children)


class DendriticTree(object):
    parameters = {}
    terminal_segments = set()
    root = None

    def __init__(self, B, E, S, N):
        """
        Parameters:
        - B: Basic branching parameter ~ expected number of branching events at an isolated segment
        - E: Size-dependency in branching ~ branching probability decrease rate coupling with number of existing segments
        - S: Order-dependency in branching ~ symmetry coefficient
        - N: total number of time bins in the full period of development
        """
        self.parameters = dict(
            B=B, E=E, S=S,
            N=N
        )
        self.root = Segment(self)
        self.terminal_segments.add(self.root)

    def grow(self, n):
        for i in range(n):
            for terminal in frozenset(self.terminal_segments):
                terminal.branch()

    
def main():
    tree = DendriticTree(B=2.5, E=0.7, S=0.5, N=50)
    tree.grow(1000)
    print("Degree at root:", tree.root.degree)
    print(tree.root.pformat())
    
    print("Function calls", counted.called)
    print("Function times", counted.timing)

if __name__ == '__main__':
    main()
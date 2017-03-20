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
    # segment's level in the dendritic tree ("centrifugal order" in [1], represented with gamma)
    order = 0
    
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

    @property
    def is_terminal(self):
        """Is this segment a terminal segment? Segment can have zero, or two
        children, being a terminal, or an intermediate segment, respectively."""
        return len(self.children) == 0

    @property
    def is_intermediate(self):
        return len(self.children) == 2

    @property
    def partition_asymmetry(self):
        """The partition asymmetry A_p at a bifurcation is defined as:
        A_p(r, s) = |r - s|/(r + s - 2), for r + r > 2
        A_p(1, 1) = 0
        """
        if self.is_intermediate:
            r, s = self.children[0].degree, self.children[1].degree
            if r + s > 2:
                return math.fabs(r - s) / (r + s - 2)
        return 0
    
    @counted
    def __init__(self, dendrite, order=0, parent=None):
        self.dendrite = dendrite
        self.parent = parent
        self.order = order

    def branching_probability_normalization_constant(self):
        s = 0
        for terminal in self.dendrite.terminal_segments:
            s += math.pow(2, -self.dendrite.parameters['S'] * terminal.order)
        return s / len(self.dendrite.terminal_segments)

    @counted
    def branch(self):
        if self.children:
            return
        
        Cinv = self.branching_probability_normalization_constant()
        B = self.dendrite.parameters['B']
        E = self.dendrite.parameters['E']
        S = self.dendrite.parameters['S']
        N = self.dendrite.parameters['N']
        n_i = len(self.dendrite.terminal_segments)
        gamma = self.order
        p_i = math.pow(2, -S * gamma) * B / Cinv / N / math.pow(n_i, E)
        
        if random.random() > p_i:
            # no branching
            return
        
        self.children = [Segment(self.dendrite, self.order + 1, self),
                         Segment(self.dendrite, self.order + 1, self)]
        self.update_degree()
        self.dendrite.terminal_segments.remove(self)
        self.dendrite.intermediate_segments.add(self)
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
        return "Segment(order=%s, degree=%s, children=[%s])" % (self.order, self.degree, children)


class DendriticTree(object):
    parameters = {}
    terminal_segments = set()
    intermediate_segments = set()
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

    @property
    def asymmetry_index(self):
        """Tree asymmetry index A_t, defined in chapter 7.1.1 in [1] is a mean
        value of all n - 1 partition asymmetries at n - 1 bifurcation points (at
        the ends of intermediate segments).
        Value is between 0 and 1, with:
            0 = complete symmetry (full binary tree)
            1 = complete asymmetry (list)
        """
        assert len(self.intermediate_segments) == len(self.terminal_segments) - 1
        if len(self.intermediate_segments) < 1:
            return 0
        s = sum([i.partition_asymmetry for i in self.intermediate_segments])
        return s / len(self.intermediate_segments)

    
def main():
    # S1-Rat Cortical Layer 2/3 Pyramidal Cell Basal Dendrites
    tree = DendriticTree(B=2.52, E=0.73, S=0.5, N=312)
    tree.grow(312)

    # Guinea Pig Purkinje Cell Dendritic Tree
    # tree = DendriticTree(B=95, E=0.69, S=-0.14, N=50)
    # tree.grow(50)

    print(tree.root.pformat())
    print("Degree at root:", tree.root.degree)
    print("Tree asymmetry index:", tree.asymmetry_index)

    print("Function calls", counted.called)
    print("Function times", counted.timing)

if __name__ == '__main__':
    main()
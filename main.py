"""
Dendritic Geometry Model, described by Jaap van Pelt et al. in [1], chapter 7.1.

[1] Computational neuroscience: Realistic modeling for experimentalists,
edited by Erik De Schutter, 2001, CRC Press.
"""
# -*- python-version >= 3.4 -*-
import math
import random
import textwrap
import statistics
from pprint import pprint

from plucky import merge

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
    
    # time of segment creation (time bin index)
    created = None

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
    def __init__(self, dendrite, order=0, parent=None, t=0):
        self.dendrite = dendrite
        self.parent = parent
        self.order = order
        self.created = t
        self.grow_initial()

    def branching_probability_normalization_constant(self):
        s = 0
        for terminal in self.dendrite.terminal_segments:
            s += math.pow(2, -self.dendrite.parameters['S'] * terminal.order)
        return s / len(self.dendrite.terminal_segments)

    @counted
    def branch(self, t):
        if self.children:
            return
        
        Cinv = self.branching_probability_normalization_constant()
        B = self.dendrite.parameters['B']
        E = self.dendrite.parameters['E']
        S = self.dendrite.parameters['S']
        N_be = self.dendrite.parameters['N_be']
        n_i = len(self.dendrite.terminal_segments)
        gamma = self.order
        p_i = math.pow(2, -S * gamma) * B / Cinv / N_be / math.pow(n_i, E)
        
        if random.random() > p_i:
            # no branching
            return
        
        self.grow_sustained(t)
        self.children = [Segment(self.dendrite, self.order + 1, self, t),
                         Segment(self.dendrite, self.order + 1, self, t)]
        self.update_degree()
        self.dendrite.terminal_segments.remove(self)
        self.dendrite.intermediate_segments.add(self)
        self.dendrite.terminal_segments.update(self.children)

    def grow_initial(self):
        """Sample from initial segment length gamma distribution, and apply to
        ``self.initial_len``."""
        p = self.dendrite.parameters
        self.initial_len = random.gammavariate(p['gamma_in'], p['beta_in']) + p['alpha_in']

    def grow_sustained(self, t):
        """Sample segment elongation rate from "branch/elongate phase" gamma distribution,
        and apply to ``self.elongated_len``."""
        p = self.dendrite.parameters
        rate = random.gammavariate(p['gamma_be'], p['beta_be']) + p['alpha_be']
        self.elongated_len = rate * (t - self.created - 1)

    def grow_only(self, dt):
        """Sample segment elongation rate from "elongate phase" gamma distribution,
        and add to ``self.elongated_len``."""
        p = self.dendrite.parameters
        rate = random.gammavariate(p['gamma_e'], p['beta_e']) + p['alpha_e']
        self.elongated_len += rate * dt

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
        return "Segment(order=%s, degree=%s, initial_len=%.4f, elongated_len=%.4f, total_len=%.4f, t_created=%d, children=[%s])" % (
            self.order, self.degree, self.initial_len, self.elongated_len, self.total_len, self.created, children
        )


class DendriticTree(object):
    parameters = None
    terminal_segments = None
    intermediate_segments = None
    root = None

    def __init__(self, B, E, S, N_be, N_e=0,
                 offset_in=0, mean_in=1, sd_in=1,
                 offset_be=0, mean_be=1, sd_be=1,
                 offset_e=0, mean_e=1, sd_e=1):
        """
        Parameters:
        - B: Basic branching parameter ~ expected number of branching events at an isolated segment
        - E: Size-dependency in branching ~ branching probability decrease rate coupling with number of existing segments
        - S: Order-dependency in branching ~ symmetry coefficient
        - N_be: total number of time bins in the branching and elongation phase
        - N_e: total number of time bins in the elongation phase
        - {offset,mean,sd}_in: gamma distribution params for initial segment lengths
        - {offset,mean,sd}_be: gamma distribution params for sustained segment elongation rate during branching and elongation phase
        - {offset,mean,sd}_e: gamma distribution params for segment elongation rate during elongation only phase
        """
        self.parameters = dict(
            B=B, E=E, S=S, N_be=N_be, N_e=N_e,

            # initial elongation, gamma distribution params:
            offset_in=offset_in, mean_in=mean_in, sd_in=sd_in,
            alpha_in=offset_in, beta_in=(sd_in**2/mean_in), gamma_in=(mean_in/sd_in)**2,

            # sustained elongation in "branching/elongation phase", gamma distribution params:
            offset_be=offset_be, mean_be=mean_be, sd_be=sd_be,
            alpha_be=offset_be, beta_be=(sd_be**2/mean_be), gamma_be=(mean_be/sd_be)**2,

            # subsequent elongation in "elongation phase", gamma distribution params:
            offset_e=offset_e, mean_e=mean_e, sd_e=sd_e,
            alpha_e=offset_e, beta_e=(sd_e**2/mean_e), gamma_e=(mean_e/sd_e)**2,
        )
        self.terminal_segments = set()
        self.intermediate_segments = set()
        self.root = Segment(self)
        self.terminal_segments.add(self.root)

    def grow(self, N_be=None, N_e=None):
        """Grow this tree: ``N_be`` iterations of branching and elongation on terminal
        segments, followed by ``N_e`` iterations of elongation only of terminal segments.
        """
        if N_be is None:
            N_be = self.parameters['N_be']
        if N_e is None:
            N_e = self.parameters['N_e']

        # branching and elongation, N_be time bins
        for i in range(N_be):
            for terminal in frozenset(self.terminal_segments):
                terminal.branch(i)

        # finish-up branching and elongation: grow un-branched terminal segments
        for terminal in self.terminal_segments:
            terminal.grow_sustained(N_be)

        # elongation only, N_e time bins
        for terminal in self.terminal_segments:
            terminal.grow_only(N_e)

    @property
    def degree(self):
        return self.root.degree

    @property
    def depth(self):
        return max([s.order for s in self.terminal_segments])

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

    @property
    def total_length(self):
        int_seg = sum([s.total_len for s in self.intermediate_segments])
        term_seg = sum([s.total_len for s in self.terminal_segments])
        return int_seg + term_seg


def simulate_and_measure(params):
    tree = DendriticTree(**params)
    tree.grow()
    return dict(
        degree=tree.degree,
        depth=tree.depth,
        asymmetry_index=tree.asymmetry_index,
        total_length=tree.total_length,
        terminal_lengths=[s.total_len for s in tree.terminal_segments],
        intermediate_lengths=[s.total_len for s in tree.intermediate_segments]
    )


def simulate(params, n):
    def add_maybe_lists(x, y):
        if not isinstance(x, list):
            x = [x]
        if not isinstance(y, list):
            y = [y]
        return x + y

    sums = {}
    for i in range(n):
        measures = simulate_and_measure(params)
        #print("Tree %d measures:" % i, measures)
        sums = merge(sums, measures, add_maybe_lists)

    stats = {}
    for k, v in sums.items():
        stats[k] = dict(total=sum(v),
                        mean=statistics.mean(v),
                        median=statistics.median(v),
                        stdev=statistics.stdev(v))
    return stats



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

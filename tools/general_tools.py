import itertools

from Bio import Seq

# DNA bases and codon->AA translation code.
bases = ['A', 'C', 'G', 'T']
codons = [''.join(cdn) for cdn in itertools.product(bases, repeat=3)]
code = {cdn: str((Seq(cdn)).translate()) for cdn in codons}
code['---'] = '-'

# Remove rare codons -- we don't want to allow these at the block junctions.
# (Phil Romero found these from various sources online.)
rare = ('ATA', 'AGG', 'TCG', 'ACG', 'AGA', 'GGA', 'CTA', 'GGG', 'CCC', 'CGG',
        'CGA')
(code.pop(r) for r in rare)

# AA->codon reverse translation code.
aminos = set(code.values())
rev_code = {a: tuple([k for k in code.keys() if code[k] == a]) for a in aminos}


class ProgressOutput:
    """ Used to show progress during iterative processes. """
    def __init__(self, total_itrs):
        """ Need total number of iterations that will be run. """
        self.total_itrs = total_itrs
        self.printed_percentage = -1

    def update(self, current_itr):
        """ Called at start of each iteration. """
        current_percentage = round(current_itr * 100 / self.total_itrs)
        if current_percentage != self.printed_percentage:
            print(f'{current_percentage}% done\r', end='', flush=True)
        self.printed_percentage = current_percentage

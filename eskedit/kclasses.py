import sys
from pyfaidx import FetchError, Fasta


class GRegion:

    def __init__(self, *args, **kwargs):
        self.default_names = ['chrom', 'start', 'stop', 'name', 'score', 'strand']
        self.fields = {'chrom': None, 'start': None, 'stop': None, 'name': None, 'score': None, 'strand': None}
        # self.fields = 6 * [None]
        for idx, key in enumerate(self.fields.keys()):
            if idx < len(args):
                self.fields[key] = str(args[idx]).strip()
        for key, value in kwargs.items():
            self.fields.update({key: value.strip()})

        try:
            self.fields['start'] = int(float(self.fields['start']))
            self.fields['stop'] = int(float(self.fields['stop']))
        except ValueError:
            raise ValueError('GRegion start and stop positions must be integers')

        self.chrom = self.fields['chrom']
        self.start = self.fields['start']
        self.stop = self.fields['stop']

    def gnomad_rep(self):
        return '{}:{}-{}'.format(self.chrom, self.start, self.stop)

    def num_fields(self):
        num_fields = 0
        for v in self.fields.values():
            if v is not None:
                num_fields += 1
        return num_fields

    def add_field(self, key_string, value):
        if not isinstance(key_string, str):
            raise KeyError('GRegion_object.add_field(key, value) must have a type \'str\' key')
        try:
            tval = self.fields.get(key_string)
            self.fields.update({key_string: value})
            return tval
        except KeyError:
            self.fields.update({key_string: value})
            return None

    def get_seq_from_fasta(self, fasta, kmer_size=1):
        # shift start left by half ksize to capture nucleotide level mutability (default is no shift)
        fa_idx_start = max(self.start - kmer_size // 2 + 1, 0)
        fa_idx_stop = self.stop + kmer_size // 2
        try:
            return fasta.get_seq(self.chrom, fa_idx_start, fa_idx_stop).seq
        except (KeyError, FetchError):
            raise KeyError('Fasta record {} : {} - {} not found'.format(self.chrom, self.start, self.stop))

    def __str__(self):
        outstring = []
        for key, value in self.fields.items():
            if value is not None:
                if key in self.default_names:
                    outstring.append(str(value))
                else:
                    outstring.append(f'{key}:{value}')
        return '\t'.join(outstring)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(other) == str(self)

    def __lt__(self, other):
        if self.fields['chrom'] == other.fields['chrom']:
            return self.fields['start'] < other.fields['start']
        else:
            # Sort numeric chromosomes numerically and non-numeric ones lexicographically
            try:
                return int(''.join(filter(str.isdigit, self.fields['chrom']))) < int(
                    ''.join(filter(str.isdigit, other.fields['chrom'])))
            except ValueError:
                return self.fields['chrom'] < other.fields['chrom']

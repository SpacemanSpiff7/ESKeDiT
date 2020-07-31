import sys
import numpy as np
from pyfaidx import FetchError, Fasta
import pandas as pd


class Model:
    def __init__(self, name: str, table: pd.DataFrame):
        self.name = name
        self.probability_dict = table.iloc[:, 0].to_dict()
        self.probability_table = table
        self.kmer_size = len(table.index[0])
        pass

    def __str__(self):
        return str({self.name: self.probability_table})

    def __repr__(self):
        return str(self)

    def get_prob(self, seq: str) -> float:
        # TODO: add to constructor or somewhere lese
        gnomad_chroms, test_chroms = 71702, 71702
        kmer_count = 0
        freq_sum = 0.0
        for start in range(len(seq) - self.kmer_size + 1):
            next_k = seq[start:start + self.kmer_size]
            if 'N' not in next_k:
                kmer_count += 1
                try:
                    freq_sum += self.probability_dict[next_k]
                except KeyError:
                    freq_sum += 0
        return freq_sum / gnomad_chroms * test_chroms


class ModelFreq:
    def __init__(self, name: str, table: pd.DataFrame):
        self.model_name = name
        self.probability_dict = table.T.to_dict()
        self.probability_table = table
        self.kmer_size = len(table.index[0])
        pass

    def __str__(self):
        return str({self.model_name: self.probability_table.head()})

    def __repr__(self):
        return str(self)

    def get_freq_array(self, kmer: str) -> np.array:
        return np.array(list(self.probability_dict.get(kmer, np.zeros(4, dtype=np.float)).values()))

    def get_expectation(self, seq: str) -> float:
        # TODO: add to constructor or somewhere lese
        gnomad_chroms, test_chroms = 71702, 71702
        scale_factor = test_chroms / gnomad_chroms
        nparray = np.array
        kmer_count = 0
        freq_sum = np.zeros(4, dtype=np.float)
        for start in range(len(seq) - self.kmer_size + 1):
            next_k = seq[start:start + self.kmer_size]
            if 'N' not in next_k:
                kmer_count += 1
                try:
                    freq_sum += nparray(list(self.probability_dict[next_k]))
                except KeyError:
                    freq_sum += 0
        return freq_sum * scale_factor


def _unpack_models(models: iter):
    new_models = []
    names = []
    for m in models:
        mname = m.model_name
        new_models.append((mname, m))
        names.append(mname)
        # self.add_model(m)
    return dict(new_models), names


class MethylationModel:
    def __init__(self, models: iter):
        self.models, self.names = _unpack_models(models)
        self.kmer_size = list(self.models.values())[0].kmer_size

    def add_model(self, model: ModelFreq):
        self.models.update({model.model_name: model})

    def get_frequency(self, model_name: str, kmer):
        return self.models.get(model_name, np.zeros(4, dtype=np.float)).get_freq_array(kmer)

    def get_methylation_frequency(self, kmer, meth_prob):
        if meth_prob < 0:
            return self.models.get('rare_transitions', np.zeros(4, dtype=np.float)).get_freq_array(kmer)
        elif meth_prob <= 0.2:
            # return low
            return self.models.get('low_methylation', np.zeros(4, dtype=np.float)).get_freq_array(kmer)
        elif 0.2 < meth_prob < 0.6:
            return self.models.get('intermediate_methylation', np.zeros(4, dtype=np.float)).get_freq_array(kmer)
        elif meth_prob > 0.6:
            return self.models.get('high_methylation', np.zeros(4, dtype=np.float)).get_freq_array(kmer)
        else:
            # throw error?
            pass


class ModelOps:
    def __init__(self):
        self.kmer_size = 0
        self.models = {}

    def add_models(self, models: iter):
        for m in models:
            if self.kmer_size == 0:
                self.kmer_size = m.kmer_size
            else:
                if self.kmer_size != m.kmer_size:
                    print(f'WARNING: K-mer size mismatch. Expected {self.kmer_size} and found {m._kmer_size}.',
                          flush=True, file=sys.stderr)
            self.add_model(m)
        return self.model_names()

    def add_model(self, model: Model):
        self.models.update({model.name: model})

    def model_names(self) -> list:
        return list(self.models.keys())

    def modeldiv(self, name_numerator: str, name_denominator: str, seq: str) -> float:
        denom = self.models[name_denominator].get_expectation(seq)
        if denom == 0:
            denom = 0.0000000000000001
        return self.models[name_numerator].get_expectation(seq) / denom

    def get_model_prob(self, name: str, seq: str) -> float:
        return self.models[name].get_expectation(seq)

    def get_model_probabilities(self, names: iter = None, seq: str = None) -> dict:
        results = {}
        if names is None:
            names = self.model_names()
        if seq is not None:
            for n in names:
                results.update({n: self.models[n].get_expectation(seq)})
        return results


# class MethylationModel(ModelOps):
#     def __init__(self, methylation_vcf: str, low=None, mid=None, hi=None, nodata=None):
#         super().__init__()

#         for name, model in {'low': low, 'mid': mid, 'hi': hi, 'nodata': nodata}.items():
#             if model is not None:
#                 self.add_model(model)


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

        if self.fields['strand'] is not None:
            self.strand = self.fields['strand']
        else:
            self.strand = None

    def gnomad_rep(self):
        return '{}:{}-{}'.format(self.chrom, self.start, self.stop)

    def strand(self):
        try:
            return self.fields['strand']
        except KeyError:
            return 'none'

    def num_fields(self):
        num_fields = 0
        for v in self.fields.values():
            if v is not None:
                num_fields += 1
        return num_fields

    def add_field(self, key_string: str, value):
        # if not isinstance(key_string, str):
        #     raise KeyError('GRegion_object.add_field(key, value) must have a type \'str\' key')
        try:
            tval = self.fields.get(key_string)
            self.fields.update({key_string: value})
            return tval
        except KeyError:
            self.fields.update({key_string: value})
            return None

    def get_seq_from_fasta(self, fasta: Fasta, kmer_size=1):
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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools as it


def generate_frequencies(transitions_path, counts_path):
    ts = pd.read_csv(transitions_path, index_col=0).sort_index()
    cts = pd.read_csv(counts_path, index_col=0).sort_index()
    cts.columns = ['counts']
    # merge to align indices
    ts = ts.merge(cts, left_on=None, right_on=None, left_index=True, right_index=True)
    ts = ts.iloc[:, :4].div(ts.iloc[:, 4], axis=0)
    return ts


def plot_individual_transitions(df):
    fig, ax = plt.subplots(4, 4)
    for r, ref in enumerate(list('ACGT')):
        for a, alt in enumerate(list('ACGT')):
            flanks = build_7mer_flanks(df, ref=ref, alt=alt)
            ax[r][a].matshow(flanks.T)


def plot_1mers(df):
    nucs = list('ACGT')

    flanks = build_1mer_flanks(df)
    level1 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 4)]
    l1_lab = nucs

    ax = plt.subplot(111)
    ax.matshow(flanks)

    ax.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax.yaxis.set_ticks_position('left')
    ax.yaxis.set_label_position('left')
    ax.set_xticks(level1)
    ax.set_xticklabels(l1_lab)
    ax.set_yticks(level1)
    ax.set_yticklabels(l1_lab)
    # ax.set_xticklabels(nucs)
    # ax.set_yticklabels(nucs)

    plt.xlabel('ALT')
    plt.ylabel('REF')

    for (i, j), z in np.ndenumerate(flanks):
        ax.text(j, i, '{:0.4f}'.format(z), ha='center', va='center')


def build_7mer_flanks(df, ref=None, alt=None, rcomp=False):
    rvcomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    nucs = list('ACGT')
    mindex = list(it.product(nucs, nucs, nucs))
    dindex = pd.MultiIndex.from_tuples(mindex)
    flanks = pd.DataFrame(index=dindex, columns=dindex)
    flanks = flanks.fillna(0)
    ct = []

    for i, r in df.iterrows():
        if ref is None:
            ref = i[3]
        if alt is None:
            alt = 0
        if i[3] == ref:
            ct.append((i, r[alt]))
        if rcomp:
            if i[3] == rvcomp[ref]:
                ct.append((i, r[rvcomp[alt]]))
    ct = np.array(ct)
    ctdf = pd.DataFrame.from_records(ct)
    ctdf.columns = ['Kmer', 'Count']
    ctdf.Count = ctdf.Count.astype('float64')
    for i, val in zip(list(ctdf.Kmer), list(ctdf.Count)):
        idx = tuple(list(i[:3]))
        col = tuple(list(i[4:]))
        flanks.loc[idx, col] = val
    return flanks.astype('float64')


def build_1mer_flanks(df, kmer_size=7):
    nucs = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    freqs = np.zeros((4, 4), dtype=np.float)
    for i, r in df.iterrows():
        from_idx = nucs.get(i[kmer_size // 2])
        freqs[from_idx] += r
    return freqs


def plot_7mer_heatmap(flanks):
    fontsize = 30
    spacing = fontsize + 5
    nucs = list('ACGT')

    plt.rcParams['figure.figsize'] = [30, 30]
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize

    level1 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 4)]
    midpt = (level1[1] - level1[0]) / 2
    level1 = [i + midpt for i in level1]
    l1_lab = nucs

    level2 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 16)]
    midpt = (level2[1] - level2[0]) / 2
    level2 = [i + midpt for i in level2]
    l2_lab = nucs * 4

    level3 = [i for i in range(0, flanks.shape[0], flanks.shape[0] // 64)]
    midpt = (level3[1] - level3[0]) / 2
    level3 = [i + midpt for i in level3]
    l3_lab = nucs * 16

    ax1 = plt.subplot(111)
    ax1.matshow(flanks.T)
    ax1.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax1.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax1.yaxis.set_ticks_position('left')
    ax1.yaxis.set_label_position('left')
    ax1.set_xticks(level3)
    ax1.set_xticklabels(l3_lab)
    ax1.set_yticks(level3)
    ax1.set_yticklabels(l3_lab)

    ax2 = ax1.twinx()  # .twinx()
    ax2.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax2.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax2.spines['bottom'].set_position(('outward', spacing))
    ax2.set_xticks(level2)
    ax2.set_xticklabels(l2_lab)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_yticks(level2)
    ax2.set_yticklabels(l2_lab)
    ax2.spines['left'].set_position(('outward', spacing))
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.set_ylim(ax1.get_ylim())

    ax3 = ax2.twinx()  # .twinx()
    ax3.xaxis.set_ticks_position('bottom')  # set the position of the second x-axis to bottom
    ax3.xaxis.set_label_position('bottom')  # set the position of the second x-axis to bottom
    ax3.spines['bottom'].set_position(('outward', spacing * 2))
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')
    ax3.spines['left'].set_position(('outward', spacing * 2))
    ax3.set_xticks(level1)
    ax3.set_xticklabels(l1_lab)
    ax3.set_yticks(level1)
    ax3.set_yticklabels(l1_lab)
    ax3.set_ylim(ax1.get_ylim())
    return ax1

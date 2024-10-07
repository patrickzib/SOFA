import sys

sys.path.insert(0, "../")

import matplotlib
import scipy as sp

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt

import warnings

warnings.simplefilter("ignore")

import numpy as np
import pandas as pd

import seaborn as sns
import scipy.fftpack as fp
from scipy.stats import norm


def readVariableLength(filepath):
    data = pd.read_csv(filepath, sep=" ", header=None)
    return np.array(data.iloc[:, 1:]), np.array(data.iloc[:, 0], dtype=np.int32)


def znorm(data):
    data = (data - np.mean(data, axis=-1, keepdims=True)) / (
            np.std(data, axis=-1, keepdims=True) + 1e-8
    )
    return data


def paa_transform(ts, n_pieces, repeat=True):
    splitted = np.array_split(ts, n_pieces)
    mean = np.asarray(list(map(lambda xs: xs.mean(axis=0), splitted)))

    n = len(ts)
    m = n / n_pieces

    if (repeat):
        means_ext = np.repeat(mean, m, axis=0)
        mean = means_ext

    return np.array(mean)


def sax_transform(ts, n_pieces, alphabet):
    alphabet_sz = len(alphabet)
    breakpoints = norm.ppf(
        np.linspace(1. / alphabet_sz, 1 - 1. / alphabet_sz, alphabet_sz - 1))

    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < breakpoints[0]
                            else (alphabet[-1] if ts_value > breakpoints[-1]
                                  else alphabet[
            np.where(breakpoints <= ts_value)[0][-1] + 1]))
                           for ts_value in ts_values])

    paa_ts = paa_transform(ts, n_pieces, False)
    # print("\tPAA", paa_ts)
    sax_ts = np.apply_along_axis(translate, 0, paa_ts), breakpoints
    # print("\tSAX", sax_ts)
    return sax_ts, paa_ts, breakpoints


def translate(ts_value, alphabet, breakpoints):
    if ts_value < breakpoints[0]:
        return 0
    elif ts_value > breakpoints[-1]:
        return len(alphabet) - 1
    else:
        return np.where(breakpoints <= ts_value)[0][-1] + 1


def dft(ts, segments):
    dft = fp.rfft(ts)

    # determine variance
    dft_variance = np.var(dft, axis=0)

    # print(dft_variance)

    # select word-length-many indices with the largest variance
    support = np.argsort(-dft_variance)[: segments + 1]

    # sort remaining indices
    support = np.sort(support)

    # support = np.arange(segments)
    return dft[:, support], support


def idft(ts, originalSize):
    return fp.irfft(ts, originalSize)


def histogram(DFTs, symbols):
    print(f"Symbols {symbols}")
    BINs = np.array([np.histogram(coeff, bins=symbols)[1] for coeff in DFTs.T])
    print(f"Symbols {BINs.shape}")

    BINs = np.zeros((DFTs.shape[1], symbols+1))
    dft = np.round(DFTs, 2)
    for letter in range(DFTs.shape[1]):
        column = np.sort(dft[:, letter])
        target_bin_width = (column[-1] - column[0]) / symbols
        BINs[letter, 0] = column[0]
        for bp in range(symbols - 1):
            BINs[letter, bp+1] = (bp + 1) * target_bin_width + column[0]

        BINs[letter, -1] = column[-1]

    return np.array(BINs)


def sfaTransform(DFTs, BINs):
    #print(DFTs.shape)
    SFAs = [np.digitize([DFTs[i]], bins=BINs[i]) for i in range(len(DFTs))]

    flat = []
    for a in SFAs:
        for b in a:
            flat.append(b)

    #print(BINs.shape)
    return np.array(flat)


def calcInverseSFA(BINs, SFA_SAMPLE, window, support):
    U_SFA_SAMPLE = SFA_SAMPLE
    U_SFA_SAMPLE[SFA_SAMPLE == len(BINs[0])] = len(BINs[0]) - 1

    upperBounds = np.array(
        [BINs[j][np.array(U_SFA_SAMPLE[j])] for j in range(len(SFA_SAMPLE))])
    lowerBounds = np.array(
        [BINs[j][np.array(SFA_SAMPLE[j] - 1)] for j in range(len(SFA_SAMPLE))])

    u_full = np.zeros((window))
    l_full = np.zeros((window))
    u_full[support] = upperBounds
    l_full[support] = lowerBounds

    upperBounds = fp.irfft(u_full)
    lowerBounds = fp.irfft(l_full)

    return (upperBounds, lowerBounds)


def plot_sax_vs_sfa(train, all_features=[4, 8, 12], ts_pos=5, alpha=256):
    alphabet = np.arange(alpha)
    # list(string.ascii_lowercase) #['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
    colors = sns.color_palette("tab10")

    print("Norming")
    train = znorm(train)
    ts1 = train[ts_pos]

    print("Done Norming")

    n = len(ts1)
    symbols = len(alphabet)

    fs_title = 32
    fs_labels = 22

    fig = plt.figure(figsize=(5 + len(all_features) * 6, 10), facecolor='white')

    ax = fig.add_subplot(2, len(all_features) + 1, 1)
    ax.set_title("Raw Time Series", fontsize=fs_title)
    ax.plot(ts1, lw=2, c=colors[0])
    ax.axis('tight')
    plt.setp(plt.xticks()[1], rotation=-45, fontsize=fs_labels)
    plt.setp(plt.yticks()[1], fontsize=fs_labels)
    ax1 = ax
    ax.set_xlabel('Time', fontsize=fs_labels)
    ax.set_ylabel('Value', fontsize=fs_labels)

    pos = 2
    for features in all_features:
        paa_ts = paa_transform(ts1, features, True)
        sax_ts, _, sax_breakpoints = sax_transform(ts1, features, alphabet)

        # append min and max value to the breakpoints
        breakpoints_new = [np.min(ts1)]
        breakpoints_new.extend(sax_breakpoints)
        breakpoints_new.append(np.max(ts1))
        breakpoints_new = np.array(breakpoints_new)

        ax = fig.add_subplot(2, len(all_features) + 1, pos, sharex=ax1, sharey=ax1)
        ax.set_title("SAX, l=" + str(features), fontsize=fs_title)
        ax.plot(ts1, lw=2, c=colors[0])
        ax.set_xlabel('Time', fontsize=fs_labels)
        ax.set_ylabel('Value', fontsize=fs_labels)

        for c, paa in enumerate(paa_ts):
            ts_value = translate(paa, alphabet, sax_breakpoints)
            maxValue = breakpoints_new[ts_value + 1]
            minValue = breakpoints_new[ts_value]
            rect = plt.Rectangle(
                (c, minValue), width=1,
                height=abs(maxValue - minValue), color=colors[1])
            ax.add_patch(rect)

        ax.plot(paa_ts, lw=2, c="black", label='PAA')

        plt.setp(plt.xticks()[1], rotation=-45, fontsize=fs_labels)
        plt.setp(plt.yticks()[1], fontsize=fs_labels)
        ax.set_xlabel('Time', fontsize=fs_labels)
        ax.set_ylabel('Value', fontsize=fs_labels)

        # text = ''.join(sax_ts[0]).upper()
        # ax.text(0, -1.6, text, fontsize=fs_title, color="black", weight='bold', family='monospace')
        plt.legend()

        print("Done SAX")

        ax = fig.add_subplot(2, len(all_features) + 1, pos + len(all_features) + 1,
                             sharex=ax1, sharey=ax1)
        ax.set_title("SFA, l=" + str(features), fontsize=fs_title)
        ax.plot(ts1, lw=2, c=colors[0], label='SFA')

        # Training
        DFTs, support = dft(train, features)
        BINs = histogram(DFTs, symbols)
        print("SFA breakpoints", np.shape(BINs))

        SFA_SAMPLE = sfaTransform(DFTs[ts_pos], BINs)
        upperBounds, lowerBounds = calcInverseSFA(BINs, SFA_SAMPLE, len(ts1), support)
        # iDFT = idft(DFTs[ts_pos], n)

        xx = np.arange(0, len(ts1))
        plt.plot(xx, upperBounds, lw=2, c=colors[1])
        plt.plot(xx, lowerBounds, lw=2, c=colors[1])
        # plt.plot(iDFT, lw=2, c="black")
        ax.fill_between(xx, lowerBounds, upperBounds, edgecolor=colors[1],
                        facecolor=colors[1])
        ax.axis('tight')
        plt.setp(plt.xticks()[1], rotation=-45, fontsize=fs_labels)
        plt.setp(plt.yticks()[1], fontsize=fs_labels)
        ax.set_xlabel('Time', fontsize=fs_labels)
        ax.set_ylabel('Value', fontsize=fs_labels)

        # text = ''.join(map(str, [chr(x + 65) for x in SFA_SAMPLE[1:]]))
        # ax.text(0, -1.6, text, fontsize=fs_title, color="black", weight='bold', family='monospace')

        print("Done SFA")
        pos = pos + 1

        # plt.legend()

    plt.tight_layout()
    return fig


train, labels = readVariableLength('Gun_Point/Gun_Point_TRAIN')

fig = plot_sax_vs_sfa(train, alpha=8)
# plt.savefig('plots/sax_vs_sfa2.pdf', format='pdf',
#                facecolor=fig.get_facecolor(), edgecolor='none')
plt.show()

file = 'queries/LenDB_queries.bin'
length = 256
data_type = np.float32

data = np.fromfile(file, dtype=data_type)

train = (data[:len(data) - len(data) % length].reshape(-1, length))
train = sp.stats.zscore(train, axis=1)
#print(train.shape)

fig = plot_sax_vs_sfa(train, all_features=[4, 8, 12, 16], ts_pos=10, alpha=128)
# plt.savefig('plots/sax_vs_sfa_sift1b.pdf', format='pdf',
#                 facecolor=fig.get_facecolor(), edgecolor='none')

plt.show()

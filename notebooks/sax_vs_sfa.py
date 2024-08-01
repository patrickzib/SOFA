# %load_ext autoreload
# %autoreload 2

import sys
sys.path.insert(0, "../")

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#%matplotlib inline
#%config InlineBackend.figure_formats = {'png', 'retina'}

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
    breakpoints = norm.ppf(np.linspace(1. / alphabet_sz, 1 - 1. / alphabet_sz, alphabet_sz - 1))

    def translate(ts_values):
        return np.asarray([(alphabet[0] if ts_value < breakpoints[0]
                            else (alphabet[-1] if ts_value > breakpoints[-1]
                                  else alphabet[
            np.where(breakpoints <= ts_value)[0][-1] + 1]))
                           for ts_value in ts_values])

    paa_ts = paa_transform(ts, n_pieces, False)
    print("\tPAA", paa_ts)
    sax_ts = np.apply_along_axis(translate, 0, paa_ts), breakpoints
    print("\tSAX", sax_ts)
    return sax_ts, paa_ts, breakpoints


def translate(ts_value, alphabet, breakpoints):
    if ts_value < breakpoints[0]:
        return 0
    elif ts_value > breakpoints[-1]:
        return len(alphabet) - 1
    else:
        return np.where(breakpoints <= ts_value)[0][-1] + 1



def dft(ts, n_pieces):
    return fp.rfft(ts)[0:n_pieces + 1]


def idft(ts, originalSize):
    return fp.irfft(ts, originalSize)


def histogram(samples, symbols):
    BINs = np.array([np.histogram(coeff, bins=symbols)[1] for coeff in DFTs.T])
    BINs[::, symbols] *= 1.01  # maxvalue
    BINs2 = []
    for i, coeff in enumerate(DFTs.T):
        BIN = []
        bin = 1
        hist = np.histogram(coeff, bins=1000)
        sums = np.sum(hist[0])
        sum = 0
        # BIN.append(hist[1][0]*5)
        for a in range(len(hist[0])):
            key = hist[0][a];
            sum = sum + key
            if sum > bin * np.ceil(sums / symbols):
                value = hist[1][a];
                BIN.append(value)
                bin = bin + 1

        BIN.append(hist[1][a])
        BINs2.append(np.array(BIN))
    return np.array(BINs2)


def sfaTransform(DFTs, BINs):
    SFAs = [np.digitize([DFTs[i]], bins=BINs[i]) for i in range(len(DFTs))]

    flat = []
    for a in SFAs:
        for b in a:
            flat.append(b)
    return np.array(flat)


def calcInverseSFA(BINs, SFA_SAMPLE, window):
    U_SFA_SAMPLE = SFA_SAMPLE
    U_SFA_SAMPLE[SFA_SAMPLE == len(BINs[0])] = len(BINs[0]) - 1

    upperBounds = [BINs[j][np.array(U_SFA_SAMPLE[j])] for j in range(len(SFA_SAMPLE))]
    lowerBounds = [BINs[j][np.array(SFA_SAMPLE[j] - 1)] for j in range(len(SFA_SAMPLE))]

    upperBounds = fp.irfft(upperBounds, window)
    lowerBounds = fp.irfft(lowerBounds, window)

    return (upperBounds, lowerBounds)

train, labels = readVariableLength('Gun_Point/Gun_Point_TRAIN')
train.shape

alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
colors = sns.color_palette("tab10")

print("Norming")
train = znorm(train)
ts_pos = 5
ts1 = train[ts_pos]
print("Done Norming")

n = len(ts1)
symbols = len(alphabet)


fs_title = 32
fs_labels = 22

fig = plt.figure(figsize=(20, 10), facecolor='white')
all_features = [4, 8, 12]

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

    text = ''.join(sax_ts[0]).upper()
    ax.text(0, -1.6, text, fontsize=fs_title, color="black", weight='bold', family='monospace')
    plt.legend()

    print("Done SAX")

    ax = fig.add_subplot(2, len(all_features) + 1, pos + len(all_features) + 1, sharex=ax1, sharey=ax1)
    ax.set_title("SFA, l=" + str(features), fontsize=fs_title)
    ax.plot(ts1, lw=2, c=colors[0], label='SFA')

    # Training
    DFTs = np.array([dft(ts, features) for ts in train])
    BINs = histogram(DFTs, symbols)
    print("SFA breakpoints", np.shape(BINs))

    SFA_SAMPLE = sfaTransform(DFTs[ts_pos], BINs)
    upperBounds, lowerBounds = calcInverseSFA(BINs, SFA_SAMPLE, len(ts1))
    # iDFT = idft(DFTs[ts_pos], n)

    xx = np.arange(0, len(ts1))
    plt.plot(xx, upperBounds, lw=2, c=colors[1])
    plt.plot(xx, lowerBounds, lw=2, c=colors[1])
    # plt.plot(iDFT, lw=2, c="black")
    ax.fill_between(xx, lowerBounds, upperBounds, edgecolor=colors[1], facecolor=colors[1])
    ax.axis('tight')
    plt.setp(plt.xticks()[1], rotation=-45, fontsize=fs_labels)
    plt.setp(plt.yticks()[1], fontsize=fs_labels)
    ax.set_xlabel('Time', fontsize=fs_labels)
    ax.set_ylabel('Value', fontsize=fs_labels)

    text = ''.join(map(str, [chr(x + 65) for x in SFA_SAMPLE[1:]]))
    ax.text(0, -1.6, text, fontsize=fs_title, color="black", weight='bold', family='monospace')

    print("Done SFA")
    pos = pos + 1

    # plt.legend()

plt.tight_layout()
plt.savefig('plots/sax_vs_sfa2.pdf', format='pdf',
            facecolor=fig.get_facecolor(), edgecolor='none')
# plt.show()
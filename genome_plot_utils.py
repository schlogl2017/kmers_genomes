import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt


def plot_base_frequency_genome(x_data, y_data, x_label, y_label):
    base_markers = {"A": "b-",
                    "C": "r-",
                    "G": "g-",
                    "T": "y-",
                    "N": "k-"}
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111)
    y_names = []
    for y in y_data:
        y_names.append(y)
        ax.plot(x_data, y_data[y], base_markers[y], label=y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    ax.legend(y_names)
    plt.grid(True)


def base_content_slide_window(sequence, name, alphabet, window, step, plot=False):
    sequence = sequence.upper()
    bases = set(alphabet.upper())
    base_freqs = defaultdict(list)
    sizes = []
    for base in bases:
        base_freqs[base] = base_freqs.get(base, [])
    for i in range(0, len(sequence) - window + 1, step):
        subseq = sequence[i:i + window]
        assert (len(subseq) == window), 'The lenght of the subsequence must have the same size of the window'
        for base in bases:
            freq = subseq.count(base) / len(subseq) * 100
            base_freqs[base].append(round(freq, 4))
        sizes.append((i + window))
    if plot:
        plot_base_frequency_genome(sizes, base_freqs, 'Genome (kb)', 'Frequencies')
        plt.title(f"Base Distribuition in {name} genome")
        plt.savefig(f"{name}_base_freq_slidewindow_plot.pdf")
    return base_freqs, sizes


def get_gc_content_window(gc_content_window, x_label, y_label, title, as_save=False):
    """Make a plot of the GC content in from a genome window
    Input:
    List of gc contents from a slide window analysis.
    Usage:
    >get_gc_content_window(gc_no,
    "Genome window", "GC content",
    "Gc contente distribuition in H. influenza genome",
    as_save=False)
    """
    plt.figure(num=None, figsize=(15, 7), dpi=100)
    plt.plot(gc_content_window, color='#1f77b4', linestyle='--', marker='.', alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label + "(%)")
    plt.title(title)
    plt.grid()
    plt.tight_layout()
    if as_save:
        plt.savefig(title)
    else:
        plt.show()


def plot_gc_skew(gc_skew, name, save=False):
    """Plots the GC skew of a sequence."""
    x = np.array(range(0, len(gc_skew)))
    plt.figure(num=None, figsize=(15, 7), dpi=100)
    plt.rcParams.update({'font.size': 12})
    plt.plot(x, gc_skew)
    plt.xlabel('Genome length (Mbp)')
    plt.ylabel('GC skew')
    plt.title(f'GC skew in the {name} genome')
    plt.tight_layout()
    if save:
        plt.savefig(f'{name}_gcSkew.pdf')
    else:
        plt.show()


def plot_gc_at_skew(sequence, k, name):
    """Plot the AT and GC skew of a sequence
    Inputs:
    sequence - string
    k - integer representing the window length
    name - Name of the species
    Output:
    plot of the At and GC skew form a sequence
    """
    seq = sequence.upper()
    x, y, z = [], [], []
    total = gc = at = 0
    window = k
    for i in range(0, len(seq), window):
        chunk = seq[i:i + window]
        g = chunk.count('G')
        c = chunk.count('C')
        a = chunk.count('A')
        t = chunk.count('T')
        if g != 0 or c != 0:
            gc = (c - g) / (g + c)
        if a != 0 or t != 0:
            at = (a - t) / (a + t)
        total += 1
        at_skew = round(at, 3)
        gc_skew = round(gc, 3)
        x.append(total)
        y.append(at_skew)
        z.append(gc_skew)
    fig, ax = plt.subplots(figsize=(14.0, 7.0), dpi=100)
    ax.plot(x,
            y,
            color='yellow',
            marker='o',
            markerfacecolor='k',
            markersize=5,
            linestyle='dashed',
            linewidth=1)
    plt.ylabel('AT ')
    plt.title(f'{name} genome AT Skew')
    plt.show()
    fig, ax = plt.subplots(figsize=(14.0, 7.0), dpi=100)
    ax.plot(x,
            z,
            color='yellow',
            marker='o',
            markerfacecolor='k',
            markersize=5,
            linestyle='dashed',
            linewidth=1)
    plt.xlabel('Window')
    plt.ylabel('GC ')
    plt.title(f'{name} genome GC Skew')
    plt.show()

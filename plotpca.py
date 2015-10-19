#!/usr/bin/env python2.7
from intensities import *
import scipy.linalg as la


def nonconstantrows(x):
    nonconst = sp.ones(x.shape[0], dtype=bool)
    for i in xrange(0, x.shape[0]):
        if all(x[i, :] == x[i, 0]):
            nonconst[i] = False
    return nonconst


def removenonconstantrows(x):
    return x[nonconstantrows(x)]


def princomp(x):
    x -= sp.mean(x, axis=0)  # Center points at origin.
    x /= sp.std(x, axis=0)  # Scale to std unit size.
    u, s, v = la.svd(x)
    ind = np.argsort(s[::-1])  # Not sorted for some reason.
    s = s[ind]
    v = v[ind]  # Recall definition of SVD.
    pca = np.dot(x, v.T)
    # cov =  np.dot(np.dot(u, s**2), u.T)
    return pca, s, v


def scatter2(x, c=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if c is not None:
        ax.scatter(x[:, 0], x[:, 1], c=c)
    else:
        ax.scatter(x[:, 0], x[:, 1])
    return ax


r0 = np.array([0, 4, 8])
r1 = np.array([1, 5, 9])
r2 = np.array([2, 6, 10])
r3 = np.array([3, 7, 11, 12])

idx = nonconstantrows(wcl.values)
pca, s, v = princomp(wcl.values[idx, :].T)
scatter2(pca)


for label, x, y in zip(wcl.columns[r0], pca[r0, 0], pca[r2, 1]):
    plt.annotate(
        label,
        xy = (x, y), xytext = (10, 40),
        textcoords = 'offset points',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

for label, x, y in zip(wcl.columns[r1], pca[r1, 0], pca[r1, 1]):
    plt.annotate(
        label,
        xy = (x, y),
        textcoords = 'offset points',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

for label, x, y in zip(wcl.columns[r2], pca[r2, 0], pca[r2, 1]):
    plt.annotate(
        label,
        xy = (x, y), xytext = (-10, -40),
        textcoords = 'offset points',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

for label, x, y in zip(wcl.columns[r3], pca[r3, 0], pca[r3, 1]):
    plt.annotate(
        label,
        xy = (x, y), xytext = (10, -40),
        textcoords = 'offset points',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

plt.show(block=False)
plt.show()

#!/usr/bin/env python2.7
from intensities import *
import scipy.linalg as la


def removeconstantrows(x):
    nonconst = sp.ones(x.shape[0], dtype=bool)
    for i in xrange(0, x.shape[0]):
        if all(x[i, :] == x[i, 0]):
            nonconst[i] = False
    return x[nonconst, :]


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


pca, s, v = princomp(removeconstantrows(wcl.values).T)
scatter2(pca)
plt.show(block=False)
plt.show()

import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools
import numpy as np


def interp_weights(xyz, uvw):
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):

    print(values.shape)
    print(vtx.shape)
    print(wts.shape)
    print(type(values))
    print(type(vtx))
    print(type(wts))
    
    print(np.take(values, vtx))
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)


m, n, d = 3500, 1, 2
# make sure no new grid point is extrapolated
bounding_cube = np.array(list(itertools.product([0, 1], repeat=d)))
xyz = np.vstack((bounding_cube, np.random.rand(m - len(bounding_cube), d)))
f = np.random.rand(m)
uvw = np.random.rand(n, d)
vtx, wts = interp_weights(xyz, uvw)

print(interpolate(f, vtx, wts))
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi


def are_collinear(vs, tol=0.01):

    if len(vs) < 2:
        return True
    v_0 = vs[0]
    for v in vs[1:]:
        if sum(abs(np.cross(v_0, v))) > tol:
            return False
    return True


def are_coplanar(vs, tol=0.01):

    if len(vs) < 3:
        return True
    v_1 = vs[0]
    for v in vs[1:]:
        if not are_collinear([v_1, v]):
            ab = np.cross(v_1, v)
            break
    else:
        return True
    for v in vs[1:]:
        if abs(np.dot(ab, v)) > tol:
            return False
    return True


def find_neighbor_ps(ps1, ps2, min_dist, max_dist):

    if min_dist >= max_dist:
        raise ValueError("min_dist should be less than max_dist!")
    dists = cdist(ps1, ps2)
    return np.array([np.where(np.logical_and(dists[i] < max_dist, dists[i] > min_dist))[0] for i in range(len(ps1))])


def find_vneighbors(points, central_points, key=1, min_dist=-float("inf"), max_dist=float("inf")):

    """
     Parameter mod can takes values 1, 2, or 3 that correspond to the
    search for domains adjacent by vertices, edges or faces.
    """

    neighbors = {i: None for i in central_points}
    vor = Voronoi(points)
    for i in central_points:
        cp = points[i]
        region = vor.regions[vor.point_region[i]]
        if -1 in region:
            raise ValueError("The domain for \"" + str(i) + "\" point is not closed!")
        local_neighbors = []
        for j in range(len(points)):
            numb_common_vertices = len(np.intersect1d(region, vor.regions[vor.point_region[j]]))
            if i != j and numb_common_vertices >= key and min_dist < np.linalg.norm(cp - points[j]) < max_dist:
                local_neighbors.append(j)
        neighbors[i] = local_neighbors
    return neighbors


def find_closest_points(p, ps, cutoff_range=(-1, float("inf")), tol=1e-8):

    dist = cdist([p], ps)[0]
    min_dist = min(dist)
    if cutoff_range[0] < min_dist < cutoff_range[1]:
        return np.where(dist < min_dist + tol)[0]
    return None


def find_equal_ps(ps1, ps2, tol=0.1):

    dists = cdist(ps1, ps2)
    inds = np.argmin(dists, axis=1)
    inds = np.array([ind if dists[i][ind] < tol else None for i, ind in enumerate(inds)])
    return inds









# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:45:34 2013

@author: localadmin
"""

import numpy as np
from matplotlib.path import Path


def sort_poly_verts(verts):
    '''Sort polygon vertices by increasing angle from the centroid

    Assumes convex polygon
    '''
    x, y = verts.T
    x_cent, y_cent = np.mean(x), np.mean(y)
    x_vect = x - x_cent
    y_vect = y - y_cent
    theta = np.arctan2(y_vect, x_vect) + np.pi / 2
    order = np.argsort(theta)
    return verts[order]


def is_poly_complex(poly_verts):
    '''Determine if a polygon is complex (any segments cross)
    '''
    segments = gather_segments(poly_verts)
    all_segments = break_segments(segments, segments)
    if len(all_segments) > len(segments):
        return True


def combine_polys(poly_old, poly_new, mode='add'):
    '''Add, subtract, or intersect two polygons
    '''
    # TODO: handle overlapping segments better
    if not mode:
        raise ValueError
    seg_old = gather_segments(poly_old)
    seg_new = gather_segments(poly_new)
    seg_new_all = break_segments(seg_new, seg_old)
    seg_old_all = break_segments(seg_old, seg_new)
    in_old = segments_in_polygon(seg_new_all, poly_old)
    in_new = segments_in_polygon(seg_old_all, poly_new)
    if mode == 'Add':
        keep_in_old = False
        keep_in_new = False
    elif mode == 'Subtract':
        keep_in_old = True
        keep_in_new = False
    else:
        keep_in_old = True
        keep_in_new = True
    seg_new = [seg for (seg, is_in_old) in zip(seg_new_all, in_old)
                    if (is_in_old == keep_in_old)]
    seg_old = [seg for (seg, is_in_new) in zip(seg_old_all, in_new)
                    if (is_in_new == keep_in_new)]
    seg = order_segments(seg_new + seg_old).tolist()
    return seg


def order_segments(segments):
    '''Piece the segments together in order, return a list of vertices
    '''
    segments = [np.array(seg).tolist() for seg in segments]
    if not segments:
        return
    verts = [segments[0][0], segments[0][1]]
    original = range(1, len(segments))
    while original:
        match = False
        for ind, segment in enumerate(segments):
            if not ind in original:
                continue
            pt1, pt2 = segment
            if np.allclose(pt1, verts[-1]):
                verts.append(pt2)
                original.remove(ind)
                match = True
            elif np.allclose(pt2, verts[-1]):
                verts.append(pt1)
                original.remove(ind)
                match = True
        if not match:
            verts.append(verts[0])
            return np.array(verts)
    return np.array(verts)


def segments_in_polygon(segments, poly_verts):
    """Go to the centroid of the segment, see if it is in the polygon
    """
    segments = np.array(segments, dtype=float)
    pt1 = segments[:, 0]
    pt2 = segments[:, 1]
    centroid = (pt1 + pt2) / 2.
    p = Path(centroid)
    return p.contains_points(poly_verts)


def gather_segments(poly_verts):
    '''Return a list of segments based on polygon vertices
    '''
    segments = []
    for ind, vert in enumerate(poly_verts):
        if ind == len(poly_verts) - 1:
            break
        segments.append((vert, poly_verts[ind + 1]))
    return segments


def get_segment_intersect(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y):
    '''Get the intersection between two line segments

    Based on response in:
    http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    '''
    s1_x = p1_x - p0_x
    s1_y = p1_y - p0_y
    s2_x = p3_x - p2_x
    s2_y = p3_y - p2_y

    s = ((-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) /
            (-s2_x * s1_y + s1_x * s2_y + 1e-7))
    t = ((s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) /
            (-s2_x * s1_y + s1_x * s2_y + 1e-7))

    keep = (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
    x = p0_x + (t * s1_x)
    y = p0_y + (t * s1_y)
    if isinstance(x, np.ndarray):
        return np.array([x[keep], y[keep]]).T
    elif keep:
        return x, y


def break_segments(poly1_segments, poly2_segments):
    '''Break up segments if they overlap any segments in the other polygon
    '''
    poly2_segments = np.array(poly2_segments, dtype=float)
    p2x = poly2_segments[:, 0, 0]
    p2y = poly2_segments[:, 0, 1]
    p3x = poly2_segments[:, 1, 0]
    p3y = poly2_segments[:, 1, 1]
    new_segs = []
    for seg1 in poly1_segments:
        inters = get_segment_intersect(seg1[0][0], seg1[0][1],
                              seg1[1][0], seg1[1][1],
                              p2x, p2y, p3x, p3y)
        if inters.size:
            if inters.size == 1:
                new_segs.append(seg1[0], inters)
                new_segs.append(inters, seg1[1])
            else:
                pt1 = seg1[0]
                diff = inters - pt1
                dist = np.hypot(diff.T[0] ** 2, diff.T[1] ** 2)
                order = np.argsort(dist)
                inters = inters[order]
                inters = [seg1[0]] + inters.tolist() + [seg1[1]]
                for ind, inter in enumerate(inters):
                    if ind < len(inters) - 1:
                        new_segs.append(np.array([inter, inters[ind + 1]]))
        else:
            new_segs.append(seg1)
    return new_segs


def create_ellipse(max_x, max_y, min_x, min_y):
    '''Create an ellipse polygon within the given extents
    '''
    a = (max_x - min_x) / 2.
    b = (max_y - min_y) / 2.
    center = [min_x + a, min_y + b]
    rad = np.arange(31) * 12 * np.pi / 180
    x = a * np.cos(rad)
    y = b * np.sin(rad)
    return np.vstack((x, y)).T + center


def create_square(max_x, max_y, min_x, min_y):
    '''Create a square polygon within the given extents
    '''
    x = [max_x, max_x, min_x, min_x, max_x]
    y = [max_y, min_y, min_y, max_y, max_y]
    return np.vstack((x, y)).T


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    poly1_verts = create_ellipse(70, 90, 15, 10)
    poly2_verts = create_ellipse(65, 100, 10, 30)

    verts = combine_polys(poly1_verts, poly2_verts, 'Add')
    verts = np.array(verts)
    #plt.plot(poly1_verts.T[0], poly1_verts.T[1])
    #plt.plot(poly2_verts.T[0], poly2_verts.T[1])
    plt.plot(verts.T[0], verts.T[1], '--')
    plt.show()

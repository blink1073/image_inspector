# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:45:34 2013

@author: localadmin
"""

import numpy as np
from matplotlib.nxutils import points_inside_poly


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
    seg_old = gather_segments(poly_old)
    seg_new = gather_segments(poly_new)
    seg_new_all = break_segments(seg_new, seg_old)
    seg_old_all = break_segments(seg_old, seg_new)
    in_old = [segment_in_polygon(pt1, pt2, poly_old)[0]
                for (pt1, pt2) in seg_new_all]
    in_new = [segment_in_polygon(pt1, pt2, poly_new)[0]
                for (pt1, pt2) in seg_old_all]
    if mode == 'sum':
        keep_in_old = False
        keep_in_new = False
    elif mode == 'diff':
        keep_in_old = True
        keep_in_new = False
    else:
        keep_in_old = True
        keep_in_new = True
    seg_new = [seg for (seg, is_in_old) in zip(seg_new_all, in_old)
                    if (is_in_old == keep_in_old)]
    seg_old = [seg for (seg, is_in_new) in zip(seg_old_all, in_new)
                    if (is_in_new == keep_in_new)]
    return order_segments(seg_new + seg_old)


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


def segment_in_polygon(pt1, pt2, poly_verts):
    """Go to the centroid of the segment, see if it is in the polygon
    """
    pt1 = np.array(pt1, dtype=float)
    pt2 = np.array(pt2, dtype=float)
    centroid = (pt1 + pt2) / 2.
    return points_inside_poly([centroid], poly_verts)


def gather_segments(poly_verts):
    '''Return a list of segments based on polygon vertices
    '''
    segments = []
    for ind, vert in enumerate(poly_verts):
        if ind == len(poly_verts) - 1:
            break
        segments.append((vert, poly_verts[ind + 1]))
    return segments


def perpendicular(a):
    '''Return the perpendicular segment to a
    '''
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b


def is_between(a, b, c):
    '''Determine if point c is between points a and b
    '''
    cross = ((c[1] - a[1]) * (b[0] - a[0]) - (c[0] - a[0]) *
                    (b[1] - a[1]))
    if abs(cross) > 1e-6:
        return False
    dot = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1]) * (b[1] - a[1])
    if dot < 0:
        return False
    sq_len_ba = ((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) *
                        (b[1] - a[1]))
    if dot > sq_len_ba:
        return False
    return True


def segment_intersect(a1, a2, b1, b2):
    '''Find the intersection between two line segments, if it exists
    '''
    a1 = np.array(a1, dtype=float)
    a2 = np.array(a2, dtype=float)
    b1 = np.array(b1, dtype=float)
    b2 = np.array(b2, dtype=float)
    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = perpendicular(da)
    denom = np.dot(dap, db)
    num = np.dot(dap, dp)
    if denom == 0:
        raise ValueError('Segments do not cross')
    inter = (num / denom) * db + b1
    if not is_between(a1, a2, inter) or not is_between(b1, b2, inter):
        raise ValueError('Segments do not cross')
    return inter.tolist()


def break_segments(poly1_segments, poly2_segments):
    '''Break up segments if they overlap any segments in the other polygon
    '''
    new_segs = []
    for seg1 in poly1_segments:
        inters = []
        for seg2 in poly2_segments:
            try:
                inter = segment_intersect(seg1[0], seg1[1], seg2[0], seg2[1])
            except ValueError:
                pass
            else:
                inters.append(inter)
        if inters:
            if len(inters) == 1:
                new_segs.append([seg1[0], inters[0]])
                new_segs.append([seg1[1], inters[0]])
            else:
                inters = np.array(inters)
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

    poly1_verts = create_square(70, 90, 15, 24)
    poly2_verts = create_square(60, 100, 20, 30)

    segment_intersect(np.array([85, 35]), np.array([25, 35]),
                      np.array([70, 70]), np.array([70, 24]))
    verts = combine_polys(poly1_verts, poly2_verts, 'diff')
    #plt.plot(poly1_verts.T[0], poly1_verts.T[1])
    #plt.plot(poly2_verts.T[0], poly2_verts.T[1])
    plt.plot(verts.T[0], verts.T[1], '--')
    plt.show()

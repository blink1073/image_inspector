# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:29:50 2013

@author: silvester
"""
import numpy as np
import scipy.ndimage as ndi
from matplotlib.path import Path


class ROI(object):

    def __init__(self, shape, source_data, geometry, source_type='image',
                 **kwargs):
        self.shape = shape.lower()
        self.source_data = source_data
        self.source_type = source_type
        self.geometry = geometry
        if 'linewidth' in kwargs:
            self.linewidth = kwargs['linewidth']
        self.data = self._get_data()
        self.extents = self._get_extents()
        if source_type == 'xy' and not source_data is None:
            self.stats = compute_stats(self.data[:, 1])
        else:
            self.stats = compute_stats(self.data)
        self.text = self._get_text()

    def _get_text(self):
        text = str(self) + '\n'
        text += '-' * 42 + '\n'
        if not self.shape == 'point':
            if isinstance(self.data, tuple):
                # TODO: handle color data
                pass
            left_col = ['mean', 'std', 'min', 'max', 'sum', 'size']
            right_col = ['median', 'ptp', '25%', '75%', '%val', '%sat']
            for (left, right) in zip(left_col, right_col):
                try:
                    left_stat = float(self.stats[left])
                    right_stat = float(self.stats[right])
                except Exception:
                    return
                line = '{0:>6} | {1:>8.3G}   || {2:>6} | {3:>8.3G}\n'
                text += line.format(left, left_stat, right, right_stat)
        return text

    def _get_data(self):
        if self.source_data is None:
            return None
        if self.shape == 'line':
            return profile_line(self.source_data, self.geometry,
                                self.linewidth)
        elif self.shape != 'point':
            return mask_data(self.source_data, self.geometry)
        else:
            return self.source_data

    def _get_extents(self):
        if self.shape in ['lasso', 'rectangle', 'ellipse']:
            path = Path(self.geometry)
            extents = path.get_extents()
            pt = extents.p0
            width = extents.width
            height = extents.height
            if self.shape == 'ellipse':
                pt = (extents.p0 + extents.p1) / 2.
            return pt, width, height
        else:
            return self.geometry

    def __str__(self):
        if self.shape in ['rectangle', 'ellipse', 'lasso']:
            pt, wid, hgt = self.extents
            x, y = pt
            text = '{0}: ({1:.4G}, {2:.4G}), {3:.4G}, {4:.4G}'
            text = text.format(self.shape.capitalize(), x, y, wid, hgt)
        elif self.shape == 'line' and not self.source_data is None:
            p1, p2 = self.geometry
            if p1[0] > p2[0]:
                p1, p2 = p2, p1
            dx = p2[0] - p1[0]
            dy = p1[1] - p2[1]
            dist = np.hypot(dx, dy)
            angle = np.arctan2(dy, dx) * 180. / np.pi
            text = 'Line: {0:.4G}, {1} deg, ({2:.4G}, {3:.4G})'
            text = text.format(dist, int(angle), p1[0], p1[1])
        elif self.shape == 'point':
            x, y = self.geometry
            if self.data:
                text = 'Point: {0:.4G} @ ({1:.4G}, {2:.4G})'
                text = text.format(self.data, x, y)
            else:
                text = 'Point: ({0:.4G}, {1:.4G})'.format(x, y)
        else:
            text = object.__str__(self)
        return text


def profile_line(img, end_points, linewidth=1):
    """Return the intensity profile of an image measured along a scan line.

    Parameters
    ----------
    img : 2d array
        The image.
    end_points: (2, 2) list
        End points ((x1, y1), (x2, y2)) of scan line.
    linewidth: int
        Width of the scan, perpendicular to the line

    Returns
    -------
    return_value : array
        The intensity profile along the scan line. The length of the profile
        is the ceil of the computed length of the scan line.
    """
    point1, point2 = end_points
    x1, y1 = point1 = np.asarray(point1, dtype=float)
    x2, y2 = point2 = np.asarray(point2, dtype=float)
    dx, dy = point2 - point1

    # Quick calculation if perfectly horizontal or vertical (remove?)
    if x1 == x2:
        pixels = img[min(y1, y2): max(y1, y2) + 1,
                     x1 - linewidth / 2:  x1 + linewidth / 2 + 1]
        intensities = pixels.mean(axis=1)
        return intensities
    elif y1 == y2:
        pixels = img[y1 - linewidth / 2:  y1 + linewidth / 2 + 1,
                     min(x1, x2): max(x1, x2) + 1]
        intensities = pixels.mean(axis=0)
        return intensities

    theta = np.arctan2(dy, dx)
    a = dy / dx
    b = y1 - a * x1
    length = np.hypot(dx, dy)

    line_x = np.linspace(min(x1, x2), max(x1, x2), np.ceil(length))
    line_y = line_x * a + b
    y_width = abs(linewidth * np.cos(theta) / 2)
    perp_ys = np.array([np.linspace(yi - y_width,
                                    yi + y_width, linewidth) for yi in line_y])
    perp_xs = - a * perp_ys + (line_x + a * line_y)[:, np.newaxis]

    perp_lines = np.array([perp_ys, perp_xs])
    pixels = ndi.map_coordinates(img, perp_lines)
    intensities = pixels.mean(axis=1)

    return intensities


def mask_data(data, verts):
    path = Path(verts)
    if isinstance(data, tuple):
        x, y = data
        pts = np.vstack((x, y)).T
        ind = np.nonzero(path.contains_points(pts))[0]
        return pts[ind]
    else:
        x, y = np.meshgrid(np.arange(data.shape[1]),
                           np.arange(data.shape[0]))
        pts = np.vstack((x.ravel(), y.ravel())).T
        ind = np.nonzero(path.contains_points(pts))[0]
        if np.rank(data) == 2:
            return data.ravel()[ind]
        else:
            new_data = [data[index].ravel()[ind] for index in range(3)]
            return new_data


def compute_stats(data):
    if data is None or data.size == 1:
        return
    real_data = data[np.isfinite(data)]
    results = dict()
    for calc in ['mean', 'std', 'min', 'max', 'sum', 'size']:
        func = getattr(np, calc)
        results[calc] = func(real_data)
    percentiles = np.percentile(real_data, [2, 25, 50, 75, 98])
    results['25%'] = percentiles[1]
    results['median'] = percentiles[2]
    results['75%'] = percentiles[3]
    results['ptp'] = percentiles[-1] - percentiles[0]
    results['%val'] = int(float(data.size) / real_data.size * 100)
    saturated = data[data == results['max']].size
    results['%sat'] = int(float(saturated) / data.size * 100)
    return results

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:29:50 2013

@author: silvester
"""
import numpy as np
from matplotlib.path import Path
from matplotlib.lines import Line2D
from base import CanvasToolBase


class ROIToolBase(CanvasToolBase):

    def __init__(self, ax, on_move=None, on_enter=None, on_release=None,
                                            useblit=True, line_props=None):
        super(ROIToolBase, self).__init__(ax, on_move=on_move,
                                        on_enter=on_enter,
                                        on_release=on_release, useblit=useblit)
        self.shape = 'none'
        props = dict(color='r', linewidth=1, solid_capstyle='butt')
        props.update(line_props if line_props is not None else {})
        self._line = Line2D([], [], visible=False, animated=True, **props)
        self._prev_line = Line2D([], [], visible=False, animated=True,
                                 **props)
        self.verts = None
        ax.add_line(self._line)
        ax.add_line(self._prev_line)
        self._artists = [self._line, self._prev_line]

    def start(self, event):
        self._busy = True

    def finalize(self):
        self._busy = False
        self.publish_roi()
        if self.verts:
            self._prev_line.set_data(zip(*self.verts))
            self._prev_line.set_visible(False)
        self.update()
        
    def ignore(self, event):
        """Return True if event should be ignored.

        This method (or a version of it) should be called at the beginning
        of any event callback.
        """
        if self.canvas.toolbar and self.canvas.toolbar.mode:
            return True
        if hasattr(event, 'button') and event.button == 3:
            return True
        return not self.active

    @property
    def source_data(self):
        if self.ax.images:
            return self.ax.images[0].get_array()
        elif self.ax.collections:
            xy = self.ax.collections[0].get_offsets()
            return xy[:, 0], xy[:, 1]
        elif self.ax.lines:
            x, y = self.ax.lines[0].get_data()
            return x, y

    def update(self):
        if not self.verts:
            return
        self._line.set_data(zip(*self.verts))
        self._line.set_visible(True)
        self.redraw()

    @property
    def roi(self):
        if self.ax.images:
            source_type = 'image'
        else:
            source_type = 'xy'
        return ROI(self.shape, self.data, self.geometry,
                   source_type=source_type)

    def publish_roi(self):
        if not self.geometry is None:
            self.process_custom_event('roi_changed', self.roi)

    @property
    def data(self):
        raise NotImplementedError


class ROI(object):

    def __init__(self, shape, data, geometry, source_type):
        self.shape = shape.lower()
        self.geometry = geometry
        self.data = data
        self.source_type = source_type
        if source_type == 'xy' and not data is None:
            self.stats = compute_stats(data[:, 1])
        elif not data is None:
            self.stats = compute_stats(data)
        self.stat_text = self._get_stat_text()
        self.handled = False

    def _get_stat_text(self):
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

    def __str__(self):
        if self.shape in ['rectangle', 'ellipse']:
            x, y, wid, hgt = self.geometry
            text = '{0}: ({1:.4G}, {2:.4G}), {3:.4G}, {4:.4G}'
            text = text.format(self.shape.capitalize(), x, y, wid, hgt)
        elif self.shape == 'lasso':
            path = Path(self.geometry)
            extents = path.get_extents()
            x, y = extents.p0
            text = 'Lasso: ({0:.4G}, {1:.4G}), {2:.4G}, {3:.4G}'
            text = text.format(x, y, extents.width, extents.height)
        elif self.shape == 'line' and not self.data is None:
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


def compute_stats(data):
    if data is None or data.size == 1:
        return
    real_data = data[np.isfinite(data)]
    if not real_data.size:
        return
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

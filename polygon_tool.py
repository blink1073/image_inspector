# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 14:45:18 2013

@author: silvester
"""
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.patches import Path
from roi import ROIToolBase
from base import ToolHandles


class PolygonToolBase(ROIToolBase):

    def __init__(self, ax, on_move=None, on_release=None, on_enter=None,
                     useblit=True, line_props=None):
        super(PolygonToolBase, self).__init__(ax, on_move=on_move,
                                on_enter=on_enter, on_release=on_release,
                                useblit=useblit)
        self.verts = []
        self._timer = self.canvas.new_timer(interval=1000)
        self._timer.add_callback(self.on_timer)
        self._timer.start()
        self._timer_count = 0

    def on_timer(self):
        if not self._line.get_visible() or self._prev_line.get_visible():
            return
        self._line.set_linestyle('--')
        self._timer_count += 1
        if self._timer_count % 2:
            self._line.set_dashes([4, 2, 6, 4])
            self._prev_line.set_dashes([4, 3, 6, 4])
        else:
            self._line.set_dashes([4, 3, 6, 4])
            self._prev_line.set_dashes([4, 2, 6, 4])
        self.redraw()

    def finalize(self):
        super(PolygonToolBase, self).finalize()

    def start(self, event):
        super(PolygonToolBase, self).start(event)
        self.verts = [(event.xdata, event.ydata)]

    @property
    def data(self):
        path = Path(self.verts)
        source_data = self.source_data
        if isinstance(source_data, tuple):
            x, y = source_data
            if x is None:
                return
            pts = np.vstack((x, y)).T
            ind = np.nonzero(path.contains_points(pts))[0]
            return pts[ind]
        else:
            x, y = np.meshgrid(np.arange(source_data.shape[1]),
                               np.arange(source_data.shape[0]))
            pts = np.vstack((x.ravel(), y.ravel())).T
            ind = np.nonzero(path.contains_points(pts))[0]
            if np.rank(source_data) == 2:
                return source_data.ravel()[ind]
            else:
                data = [source_data[index].ravel()[ind] for index in range(3)]
                return data


class LassoSelection(PolygonToolBase):

    def __init__(self, ax, on_move=None, on_release=None, on_enter=None,
                     useblit=True, line_props=None):
        super(LassoSelection, self).__init__(ax, on_move=on_move,
                                            on_release=on_release,
                                                on_enter=on_enter,
                                                useblit=useblit,
                                                line_props=line_props)
        self.shape = 'lasso'
        self.mode = 'lasso'
        self._prev_angle = 1e6
        self._indicator = Rectangle((0, 0), 0, 0)
        self._indicator.set_visible(False)
        ax.add_patch(self._indicator)
        self._artists.append(self._indicator)

    def on_move(self, event):
        if self.verts is None or not self._busy:
            return
        if self.mode == 'polygon':
            self.verts[-1] = (event.xdata, event.ydata)
        else:
            # see if we are moving in the same direction
            point2 = (event.xdata, event.ydata)
            point1 = np.asarray(self.verts[-1], dtype=float)
            point2 = np.asarray(point2, dtype=float)
            dx, dy = point2 - point1
            theta = np.arctan2(dy, dx)
            if theta == self._prev_angle:
                self.verts[-1] = (event.xdata, event.ydata)
            else:
                self.verts.append((event.xdata, event.ydata))
            self._prev_angle = theta
        try:
            self._show_indicator()
        except TypeError:
            pass
        self.update()

    def on_mouse_release(self, event):
        self.mode = 'polygon'

    def on_mouse_press(self, event):
        if self._indicator.get_visible() or event.dblclick:
            self.verts = self.verts[:-1]
            self.finalize()
            return
        if not self._busy:
            self.start(event)
        self.mode = 'lasso'
        self.verts.append((event.xdata, event.ydata))
        self.update()

    def start(self, event):
        super(LassoSelection, self).start(event)
        self._prev_angle = 1e6

    def finalize(self):
        self.verts.append(self.verts[0])
        self._indicator.set_visible(False)
        super(LassoSelection, self).finalize()

    def _show_indicator(self):
        if not self.verts:
            return
        bounds = self.ax.dataLim.bounds
        wid = float(bounds[2] - bounds[0])
        hgt = float(bounds[3] - bounds[1])
        orig = self.verts[0]
        curr = self.verts[-1]
        # see if we are within 2% in x and y
        if (abs(curr[0] - orig[0]) / wid < 0.02 and
            abs(curr[1] - orig[1]) / hgt < 0.02):
                wid /= 100.
                hgt /= 100.
                self._indicator.set_xy((orig[0] - wid / 2.,
                                              orig[1] - hgt / 2.))
                self._indicator.set_width(wid * 2)
                self._indicator.set_height(hgt * 2)
                if not self._indicator.get_visible():
                    self._indicator.set_visible(True)
                    self.redraw()
        elif self._indicator.get_visible():
            self._indicator.set_visible(False)
            self.redraw()

    @property
    def geometry(self):
        return self.verts

    @geometry.setter
    def geometry(self, verts):
        self.verts = verts
        self.update()


class RectangleSelection(PolygonToolBase):

    def __init__(self, ax, maxdist=10, on_move=None,
                 on_release=None, on_enter=None, useblit=True,
                 line_props=None):
        super(RectangleSelection, self).__init__(ax, on_move=on_move,
                                            on_release=on_release,
                                                on_enter=on_enter,
                                                useblit=useblit,
                                                line_props=line_props)
        self.anchor = None
        self.origin = None
        self.origin_pix = None
        self.shape = 'rectangle'
        self.maxdist = maxdist
        self._extents = [0, 0, 0, 0]
        self.active_handle = None
        self.modifiers = set()
        x = (0, 0, 0, 0)
        y = (0, 0, 0, 0)
        props = dict()
        self._center_handle = ToolHandles(ax, x, y, marker='s',
                                          marker_props=props)
        self._corner_handles = ToolHandles(ax, x, y, marker_props=props)
        self._edge_handles = ToolHandles(ax, x, y, marker='s',
                                         marker_props=props)
        self._corner_order = ['NW', 'NE', 'SE', 'SW']
        self._edge_order = ['W', 'N', 'E', 'S']
        self._artists.extend([self._center_handle.artist,
                         self._corner_handles.artist,
                         self._edge_handles.artist])
        self.connect_event('key_press_event', self.on_key)
        self.connect_event('key_release_event', self.off_key)

    def start(self, event):
        super(RectangleSelection, self).start(event)
        self.anchor = None
        self.origin = self.verts[0]
        self.origin_pix = (event.x, event.y)
        if not self.active_handle == 'C':
            self.extents = [0, 0, 0, 0]
        self._line.set_visible(True)

    def on_key(self, event):
        if self._busy:
            self.modifiers.add(event.key)

    def off_key(self, event):
        if self._busy and event.key in self.modifiers:
            self.modifiers.remove(event.key)

    def on_move(self, event):
        if not self._busy or not event.xdata:
            return
        if self.active_handle and not self.active_handle == 'C':
            x1, x2, y1, y2 = self._extents_on_press
            if self.active_handle in ['E', 'W'] + self._corner_order:
                x2 = event.xdata
            if self.active_handle in ['N', 'S'] + self._corner_order:
                y2 = event.ydata
            self.modifiers = set()
            self.extents = np.array([x1, x2, y1, y2])
        elif ' ' in self.modifiers or self.active_handle == 'C':
            # move command
            if not self.anchor:
                self.anchor = event.xdata, event.ydata
            else:
                dx = event.xdata - self.anchor[0]
                dy = event.ydata - self.anchor[1]
                x1, x2, y1, y2 = self.extents
                x1 += dx
                x2 += dx
                y1 += dy
                y2 += dy
                self.extents = np.array([x1, x2, y1, y2])
                self.origin = [self.origin[0] + dx, self.origin[1] + dy]
                self.anchor = event.xdata, event.ydata
        else:
            dx = (event.xdata - self.origin[0]) / 2.
            dy = (event.ydata - self.origin[1]) / 2.
            center = np.array(self.origin)
            if 'shift' in self.modifiers or 'ctrl+shift' in self.modifiers:
                # square command
                dx_pix = abs(event.x - self.origin_pix[0])
                dy_pix = abs(event.y - self.origin_pix[1])
                if not dx_pix:
                    return
                maxd = max(abs(dx_pix), abs(dy_pix))
                if abs(dx_pix) < maxd:
                    dx *= maxd / abs(dx_pix)
                if abs(dy_pix) < maxd:
                    dy *= maxd / abs(dy_pix)
            if 'control' in self.modifiers or 'ctrl+shift' in self.modifiers:
                # from center
                dx *= 2
                dy *= 2
            else:
                # from corner
                center += np.array([dx, dy])
            self.extents = np.array([center[0] - dx, center[0] + dx,
                                    center[1] - dy, center[1] + dy])

    def set_verts(self):
        x1, x2, y1, y2 = self.extents
        self.verts = [[x1, y1], [x1, y2], [x2, y2], [x2, y1],
                          [x1, y1]]
        self.update()

    def on_mouse_release(self, event):
        self._extents_on_press = None
        self.modifiers = set()
        self.finalize()

    def on_mouse_press(self, event):
        self._set_active_handle(event)
        self.start(event)

    def _set_active_handle(self, event):
        """Set active handle based on the location of the mouse event"""
        # Note: event.xdata/ydata in data coordinates, event.x/y in pixels
        c_idx, c_dist = self._corner_handles.closest(event.x, event.y)
        e_idx, e_dist = self._edge_handles.closest(event.x, event.y)
        m_idx, m_dist = self._center_handle.closest(event.x, event.y)

        # Set active handle as closest handle, if mouse click is close enough.
        if c_dist > self.maxdist and e_dist > self.maxdist:
            self.active_handle = None
        elif c_dist < e_dist:
            self.active_handle = self._corner_order[c_idx]
        else:
            self.active_handle = self._edge_order[e_idx]
        if not self.active_handle and m_dist < self.maxdist:
            self.active_handle = 'C'
        # Save coordinates of rectangle at the start of handle movement.
        x1, x2, y1, y2 = self.extents
        # Switch variables so that only x2 and/or y2 are updated on move.
        if self.active_handle in ['W', 'SW', 'NW']:
            x1, x2 = x2, event.xdata
        if self.active_handle in ['N', 'NW', 'NE']:
            y1, y2 = y2, event.ydata
        self._extents_on_press = x1, x2, y1, y2

    @property
    def _rect_bbox(self):
        if not len(self.extents):
            return 0, 0, 0, 0
        x1, x2, y1, y2 = self.extents.tolist()
        x0, x1 = np.sort((x1, x2))
        y0, y1 = np.sort((y1, y2))
        width = x1 - x0
        height = y1 - y0
        return x0, y0, width, height

    @property
    def corners(self):
        """Corners of rectangle from lower left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        xc = x0, x0 + width, x0 + width, x0
        yc = y0, y0, y0 + height, y0 + height
        return xc, yc

    @property
    def edge_centers(self):
        """Midpoint of rectangle edges from left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        w = width / 2.
        h = height / 2.
        xe = x0, x0 + w, x0 + width, x0 + w
        ye = y0 + h, y0, y0 + h, y0 + height
        return xe, ye

    @property
    def extents(self):
        return self._extents

    @extents.setter
    def extents(self, extents):
        x1, x2, y1, y2 = extents
        xmin, xmax = np.sort([x1, x2])
        ymin, ymax = np.sort([y1, y2])
        self._extents = np.array([xmin, xmax, ymin, ymax])

        # Update displayed handles
        self._center_handle.set_data((x2 + x1) / 2, (y2 + y1) / 2)
        self._corner_handles.set_data(*self.corners)
        self._edge_handles.set_data(*self.edge_centers)

        self._center_handle.set_visible(True)
        self._corner_handles.set_visible(True)
        self._edge_handles.set_visible(True)
        self.set_verts()

    @property
    def geometry(self):
        return self._rect_bbox

    @geometry.setter
    def geometry(self, x0, y0, width, height):
        self.extents = [x0, x0 + width, y0, y0 + height]


class EllipseSelection(RectangleSelection):

    def __init__(self, ax, maxdist=10, on_move=None,
                 on_release=None, on_enter=None, useblit=True,
                 line_props=None):
        super(EllipseSelection, self).__init__(ax, on_move=on_move,
                                            on_release=on_release,
                                                on_enter=on_enter,
                                                maxdist=maxdist, useblit=True,
                                                line_props=line_props)
        self.shape = 'ellipse'

    def set_verts(self):
        x1, x2, y1, y2 = self.extents
        # create an elliptical path within the extents
        min_x, max_x = np.sort((x1, x2))
        min_y, max_y = np.sort((y1, y2))
        a = (max_x - min_x) / 2.
        b = (max_y - min_y) / 2.
        center = [min_x + a, min_y + b]
        rad = np.arange(31) * 12 * np.pi / 180
        x = a * np.cos(rad)
        y = b * np.sin(rad)
        self.verts = np.vstack((x, y)).T + center
        self.verts = self.verts.tolist()
        self.update()

    @property
    def geometry(self):
        x0, y0, width, height = self._rect_bbox
        return x0 + width / 2., y0 + width / 2., width, height

    @geometry.setter
    def geometry(self, xc, yc, width, height):
        self.extents = [xc - width / 2., yc - height / 2., width, height]


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from skimage import data

    image = data.camera()

    f, ax = plt.subplots()
    ax.imshow(image, interpolation='nearest')
    h, w = image.shape

    def roi_changed(roi):
        print roi.stat_text

    tool = RectangleSelection(ax)
    tool.connect_event('roi_changed', roi_changed)
    plt.show()

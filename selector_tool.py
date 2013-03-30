# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 21:11:25 2013

@author: silvester
"""
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import polygon_math

from base import CanvasToolBase, ToolHandles
from roi import ROI


class SelectionTool(CanvasToolBase):
    """Widget for selecting a rectangular region in a plot.

    After making the desired selection, press "Enter" to accept the selection
    and call the `on_enter` callback function.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Matplotlib axes where tool is displayed.
    on_move : function
        Function called whenever a control handle is moved.
        This function must accept the rectangle extents as the only argument.
    on_finish : function
        Function called whenever the control handle is released.
    on_enter : function
        Function called whenever the "enter" key is pressed.
    maxdist : float
        Maximum pixel distance allowed when selecting control handle.
    rect_props : dict
        Properties for :class:`matplotlib.patches.Rectangle`. This class
        redefines defaults in :class:`matplotlib.widgets.RectangleSelector`.

    Attributes
    ----------
    vertices : list of tuples
        (x, y) points definining the shape.
    """
    def __init__(self, ax, on_move=None, on_release=None, on_enter=None,
                 on_finish=None, maxdist=10, lineprops=None, shape='rectangle'):
        CanvasToolBase.__init__(self, ax, on_move=on_move,
                                on_enter=on_enter, on_release=on_finish)

        if on_enter is None:
            def on_enter(vertices):
                print(vertices)
        self.callback_on_enter = on_enter

        props = dict(linestyle='-.', color='black')
        props.update(lineprops if lineprops is not None else {})

        self._line = Line2D([], [], **props)
        self._prev_line = Line2D([], [], **props)
        for line in [self._line, self._prev_line]:
            line.set_visible(False)
            self.ax.add_line(line)

        self.connect_event('button_press_event', self.onpress)
        self.connect_event('button_release_event', self.onrelease)
        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('key_press_event', self.onkey)
        self.connect_event('key_release_event', self.offkey)

        self.mode = 'New'
        self.callback_on_finish = self.callback_on_release

        self._prev_verts = None
        self._prev_prev_verts = None

        self.modifiers = set()
        self._prev_modifiers = set()
        self.drawing = False

        self._timer = self.ax.figure.canvas.new_timer(interval=1000)
        self._timer.add_callback(self.ontimer)
        self._timer.start()
        self._timer_count = 0

        self._lasso_tool = LassoSelection(ax)
        self._marquee_tool = MarqeeSelection(ax, maxdist=maxdist)

        self.set_shape(shape)

        self._artists = [self._line,
                         self._prev_line]
        self._artists.extend(self._lasso_tool._artists)
        self._artists.extend(self._marquee_tool._artists)

    def onkey(self, event):
        '''Update our modifiers on a key press
        '''
        if self.tool.active:
            self.tool.add_modifier(event.key)
        else:
            self.modifiers.add(event.key)
        if event.key == 'ctrl+r':
            self.set_shape('Rectangle')
        elif event.key == 'ctrl+l':
            self.set_shape('Lasso')
        elif event.key == 'ctrl+e':
            self.set_shape('Ellipse')

    def offkey(self, event):
        '''Update our modifiers on a key release
        '''
        if event.key in self.tool.modifiers:
            self.tool.remove_modifier(event.key)
        elif not self.tool.active and event.key in self.modifiers:
            self.modifiers.remove(event.key)

    def onmove(self, event):
        '''Update our verts on a mouse movement
        '''
        if self.ignore(event):
            return
        self.tool.onmove(event)
        self.update(self.tool.verts)

    def onrelease(self, event):
        '''React to a release of the mouse button
        '''
        if self.ignore(event):
            return
        self.tool.onrelease(event)
        self.update(self.tool.verts)

    def onpress(self, event):
        '''React to the pressing of a mouse button
        '''
        if self.ignore(event):
            return
        self.drawing = True
        if self.modifiers:
            self._prev_line.set_visible(True)
        self.tool.onpress(event)
        if hasattr(self.tool, 'active_handle') and self.tool.active_handle:
            self.modifiers = self.modifiers.union(self._prev_modifiers)
            self._prev_verts = self._prev_prev_verts
        self.update(self.tool.verts)

    def ignore(self, event):
        if not hasattr(event, 'inaxes'):
            return False
        if self.active and (event.inaxes == self.ax):
            return False
        else:
            return True

    def update(self, event):
        if self.tool.finished and not self.tool.verts is None:
            self.finalize()
            return
        self.draw_line(self.tool.verts)

    def draw_line(self, verts):
        if not verts:
            return
        self._line.set_data(zip(*verts))
        self._line.set_visible(True)
        if len(verts) == 1:
            self.tool.set_visible(False)
        self.redraw()

    def ontimer(self):
        if not self._line.get_visible():
            return
        self._line.set_linestyle('--')
        self._timer_count += 1
        if self._timer_count % 2:
            self._line.set_dashes([4, 2, 6, 4])
        else:
            self._line.set_dashes([4, 3, 6, 4])
        self.redraw()

    def finalize(self):
        '''Take the appropriate action based on mode
        '''
        if not self.drawing:
            return
        verts = self.tool.verts
        if self._prev_verts is not None and self.modifiers:
            if 'ctrl+shift' in self.modifiers:
                mode = 'Intersect'
            elif 'shift' in self.modifiers:
                mode = 'Add'
            elif 'control' in self.modifiers:
                mode = 'Subtract'
            else:
                mode = None
            try:
                # convert to axis points
                old = self.ax.transAxes.transform(self._prev_verts)
                new = self.ax.transAxes.transform(self.tool.verts)
                verts = polygon_math.combine_polys(old, new, mode)
                verts = self.ax.transAxes.inverted().transform(verts).tolist()
            except (ValueError, IndexError):
                pass
        self.callback_on_finish(verts)
        self._prev_prev_verts = self._prev_verts
        self._prev_verts = verts
        self.tool.active = False

        self.draw_line(verts)
        self._prev_line.set_data(zip(*verts))
        self._prev_line.set_visible(False)
        self._prev_modifiers = self.modifiers
        self.modifiers = set()
        self.drawing = False
        self.tool.cleanup()

    def set_shape(self, shape):
        if shape.lower() == 'rectangle':
            self.tool = self._marquee_tool
            self.tool.shape = 'Rectangle'
        elif shape.lower() == 'ellipse':
            self.tool = self._marquee_tool
            self.tool.shape = 'Ellipse'
        elif shape.lower() == 'lasso':
            self.tool = self._lasso_tool
        self._marquee_tool.set_visible(False)
        self.shape = shape

    @property
    def geometry(self):
        return self.tool.verts

    @geometry.setter
    def geometry(self, shape, points):
        self.set_shape(shape)
        self.tool.verts = points
        self.finalize()

    def activate(self):
        self.active = True
        if hasattr(self, 'tool'):
            self.tool.active = True
        self.redraw()
        self.modifiers = set()

    def deactivate(self):
        self.active = False
        self.tool.active = False
        self.set_visible(False)
        self.modifiers = set()

    @property
    def data(self):
        return self.tool.data

    @property
    def roi(self):
        return self.tool.roi


class BaseSelector(CanvasToolBase):

    def __init__(self, ax):
        CanvasToolBase.__init__(self, ax)
        self.verts = None
        self.modifiers = set()
        self.finished = False
        self._prev_verts = None
        self.shape = 'rectangle'

    def onmove(self, event):
        pass

    def finalize(self):
        self.finished = True
        self.modifiers = set()
        self.active = False
        if not self.verts is None and len(self.verts) > 1:
            self.canvas.callbacks.process('roi_changed', self.roi)

    def onrelease(self, event):
        self.finalize()

    def onpress(self, event):
        self.start(event)

    def add_modifier(self, modifier):
        self.modifiers.add(modifier)

    def remove_modifier(self, modifier):
        self.modifiers.remove(modifier)

    def start(self, event):
        self.active = True
        self.verts = [(event.xdata, event.ydata)]
        self.finished = False

    def cleanup(self):
        pass

    @property
    def geometry(self):
        return self.verts

    @property
    def data(self):
        if self.ax.images:
            data = self.ax.images[0].get_array()
            return data
        elif self.ax.lines:
            x, y = self.ax.lines[0].get_data()
            return x, y

    @property
    def roi(self):
        return ROI(self.shape, self.data, self.geometry)


class LassoSelection(BaseSelector):

    def __init__(self, ax):
        BaseSelector.__init__(self, ax)
        self.shape = 'lasso'
        self.mode = 'lasso'
        self._indicator = Rectangle((0, 0), 0, 0)
        self._indicator.set_visible(False)
        self._prev_angle = 1e6
        ax.add_patch(self._indicator)
        self._artists = [self._indicator]

    def onmove(self, event):
        if not self.active or self.verts is None:
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
        self.show_close()

    def onrelease(self, event):
        if not self.active:
            return
        self.mode = 'polygon'

    def onpress(self, event):
        if self._indicator.get_visible() or event.dblclick:
            self.verts = self.verts[:-1]
            self.finalize()
            return
        if not self.verts:
            self.start(event)
        self.mode = 'lasso'
        self.verts.append((event.xdata, event.ydata))

    def finalize(self):
        self.verts.append(self.verts[0])
        self._prev_angle = 1e6
        BaseSelector.finalize(self)

    def cleanup(self):
        self.verts = None
        self._indicator.set_visible(False)

    def show_close(self):
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


class MarqeeSelection(BaseSelector):

    def __init__(self, ax, shape='Rectangle', maxdist=10):
        BaseSelector.__init__(self, ax)
        self.anchor = None
        self.origin = None
        self.origin_pix = None
        self.extents = [0, 0, 0, 0]
        self.shape = shape
        self.maxdist = maxdist
        self.active_handle = None
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
        self._artists = [self._center_handle.artist,
                         self._corner_handles.artist,
                         self._edge_handles.artist]

    def start(self, event):
        self.verts = None
        self.anchor = None
        BaseSelector.start(self, event)
        self.origin = self.verts[0]
        self.origin_pix = (event.x, event.y)
        if not self.active_handle == 'C':
            self.extents = [0, 0, 0, 0]
        self.set_visible(True)

    def onmove(self, event):
        if not self.origin or not self.active or self.finished:
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
            if 'shift' in self.modifiers:
                # square command
                dx_pix = abs(event.x - self.origin_pix[0])
                dy_pix = abs(event.y - self.origin_pix[1])
                maxd = max(abs(dx_pix), abs(dy_pix))
                if abs(dx_pix) < maxd:
                    dx *= maxd / abs(dx_pix)
                if abs(dy_pix) < max:
                    dy *= maxd / abs(dy_pix)
            if 'control' in self.modifiers:
                # from center
                dx *= 2
                dy *= 2
            else:
                # from corner
                center += np.array([dx, dy])
            self.extents = np.array([center[0] - dx, center[0] + dx,
                                    center[1] - dy, center[1] + dy])
        self.set_verts()

    def set_verts(self):
        self.set_extents()
        x1, x2, y1, y2 = self.extents
        if self.shape.lower() == 'rectangle':
            self.verts = [[x1, y1], [x1, y2], [x2, y2], [x2, y1],
                          [x1, y1]]
        else:
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

    def onrelease(self, event):
        self._extents_on_press = None
        self.finalize()

    def onpress(self, event):
        self._set_active_handle(event)
        BaseSelector.onpress(self, event)

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

    def set_extents(self):
        x1, x2, y1, y2 = self.extents
        xmin, xmax = np.sort([x1, x2])
        ymin, ymax = np.sort([y1, y2])
        self.extents = np.array([xmin, xmax, ymin, ymax])

        # Update displayed handles
        self._center_handle.set_data((x2 + x1) / 2, (y2 + y1) / 2)
        self._corner_handles.set_data(*self.corners)
        self._edge_handles.set_data(*self.edge_centers)

        self._center_handle.set_visible(True)
        self._corner_handles.set_visible(True)
        self._edge_handles.set_visible(True)
        self.redraw()


if __name__ == '__main__':
    np.testing.rundocs()
    import matplotlib.pyplot as plt
    from skimage import data

    image = data.camera()

    def print_geometry(geometry):
        print geometry

    def roi_changed(roi):
        print roi.shape, len(roi.geometry), type(roi.data)

    f, ax = plt.subplots()
    ax.imshow(image, interpolation='nearest')
    #ax.plot(np.random.random(10))
    tool = SelectionTool(ax, shape='lasso', on_release=print_geometry)
    tool.connect_event('roi_changed', roi_changed)
    plt.show()

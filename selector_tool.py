# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 21:11:25 2013

@author: silvester
"""
import numpy as np
from matplotlib.widgets import AxesWidget
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.path import Path
import polygon_math


class Transform(object):
    pass


# TODO: finishing a current tool brings out a transform tool for rectangle and square
#        but it only works on the current tool


class Selector(AxesWidget):

    def __init__(self, ax, onselect=None, useblit=True, lineprops=None,
                 shape='Lasso'):
        AxesWidget.__init__(self, ax)

        self.useblit = useblit
        self.onselect = onselect

        if lineprops is None:
            lineprops = dict()
        lineprops.setdefault('linestyle', '-.')
        lineprops['color'] = 'black'
        self.line = Line2D([], [], **lineprops)
        self.line.set_visible(False)
        self.ax.add_line(self.line)

        self.prev_line = Line2D([], [], **lineprops)
        self.prev_line.set_visible(False)
        self.ax.add_line(self.prev_line)

        lineprops['linestyle'] = '-'
        self.extent_line = Line2D([], [], **lineprops)
        self.extent_line.set_visible(False)
        self.ax.add_line(self.extent_line)

        self.connect_event('button_press_event', self.onpress)
        self.connect_event('button_release_event', self.onrelease)
        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.update_background)
        self.connect_event('key_press_event', self.onkey)
        self.connect_event('key_release_event', self.offkey)

        self.mode = 'New'
        self.set_shape(shape)
        self.ax.selection = []

        self._patch = Rectangle((0, 0), 0, 0)
        self._patch.set_visible(False)
        self.ax.add_patch(self._patch)
        self.modifiers = None

        self.timer = self.ax.figure.canvas.new_timer(interval=1000)
        self.timer.add_callback(self.ontimer)
        self.timer.start()
        self.timer_count = 0

    def onkey(self, event):
        '''Update our modifiers on a key press
        '''
        if self.tool.active:
            self.tool.add_modifier(event.key)
        else:
            self.modifiers.add(event.key)

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
        self.tool.onpress(event)
        self.prev_line.set_visible(True)
        if (isinstance(self.tool, LassoSelection) and
                (self._patch.get_visible() or event.dblclick)):
            self.tool.verts = self.tool.verts[:-2]
            self.tool.finalize()
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
        self.show_close()
        self.draw_line(self.tool.verts)

    def draw_line(self, verts):
        if not verts:
            return
        self.line.set_data(zip(*verts))
        self.line.set_visible(True)

        if isinstance(self.tool, EllipticalSelection):
            self.extent_line.set_visible(True)
            ext_verts = self.tool.extents.tolist() + [self.tool.extents[0]]
            self.extent_line.set_data(zip(*ext_verts))
        self.update_plot()

    def update_plot(self):
        if self.useblit:
            self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.line)
            self.ax.draw_artist(self.prev_line)
            self.ax.draw_artist(self.extent_line)
            self.ax.draw_artist(self._patch)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()

    def ontimer(self):
        if not self.line.get_visible():
            return
        self.line.set_linestyle('--')
        self.timer_count += 1
        if self.timer_count % 2:
            self.line.set_dashes([4, 2, 6, 4])
        else:
            self.line.set_dashes([4, 3, 6, 4])
        self.update_plot()

    def finalize(self):
        '''Take the appropriate action based on mode
        '''
        verts = self.tool.verts
        if self.ax.selection and self.modifiers:
            if 'ctrl+shift' in self.modifiers:
                mode = 'Intersect'
            elif 'shift' in self.modifiers:
                mode = 'Add'
            elif 'control' in self.modifiers:
                mode = 'Subtract'
            else:
                mode = None
            try:
                verts = polygon_math.combine_polys(self.ax.selection,
                                                   self.tool.verts,
                                                   mode)
            except (ValueError, IndexError):
                pass
        self.onselect(verts)
        self.ax.selection = verts
        self.tool.active = False
        self.tool.verts = None
        self._patch.set_visible(False)
        self.draw_line(verts)
        self.prev_line.set_data(zip(*verts))
        self.prev_line.set_visible(False)
        self.modifiers = set()

    def set_shape(self, shape):
        if shape == 'Rectangle':
            self.tool = RectangleSelection()
        elif shape == 'Ellipse':
            self.tool = EllipticalSelection()
        elif shape == 'Lasso':
            self.tool = LassoSelection()

    def update_background(self, event):
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def show_close(self):
        if not isinstance(self.tool, LassoSelection) or not self.tool.verts:
            return
        bounds = self.ax.dataLim.bounds
        wid = float(bounds[2] - bounds[0])
        hgt = float(bounds[3] - bounds[1])
        orig = self.tool.verts[0]
        curr = self.tool.verts[-1]
        # see if we are within 2% in x and y
        if (abs(curr[0] - orig[0]) / wid < 0.02 and
            abs(curr[1] - orig[1]) / hgt < 0.02):
                wid /= 100.
                hgt /= 100.
                self._patch.set_xy((orig[0] - wid / 2., orig[1] - hgt / 2.))
                self._patch.set_width(wid)
                self._patch.set_height(hgt)
                self._patch.set_visible(True)
                self.update_plot()
        else:
            self._patch.set_visible(False)
            self.update_plot()


class SelectorTool(object):

    def __init__(self):
        self.active = False
        self.verts = []
        self.modifiers = set()
        self.finished = False

    def onmove(self, event):
        pass

    def finalize(self):
        self.finished = True

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


class LassoSelection(SelectorTool):

    def __init__(self):
        SelectorTool.__init__(self)
        self.mode = 'lasso'

    def onmove(self, event):
        if not self.active or not len(self.verts):
            return
        if self.mode == 'polygon':
            self.verts[-1] = (event.xdata, event.ydata)
        else:
            self.verts.append((event.xdata, event.ydata))

    def onrelease(self, event):
        if not self.active:
            return
        self.mode = 'polygon'

    def onpress(self, event):
        if not self.verts:
            self.start(event)
        self.mode = 'lasso'
        self.verts.append((event.xdata, event.ydata))

    def finalize(self):
        self.verts.append(self.verts[0])
        self.finished = True


class RectangleSelection(SelectorTool):

    def __init__(self):
        SelectorTool.__init__(self)
        self.anchor = None
        self.origin = None
        self.origin_pix = None
        self.extents = []

    def start(self, event):
        SelectorTool.start(self, event)
        self.origin = self.verts[0]
        self.origin_pix = (event.x, event.y)
        self.extents = np.vstack(self.verts * 4)

    def onmove(self, event):
        if not self.origin or not self.active:
            return
        if ' ' in self.modifiers:
            # move command
            if not self.anchor:
                self.anchor = event.xdata, event.ydata
            else:
                dx = event.xdata - self.anchor[0]
                dy = event.ydata - self.anchor[1]
                self.extents += np.array([dx, dy])
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
            self.extents = np.vstack((center + [dx, dy],
                                       center + [dx, -dy],
                                       center + [-dx, -dy],
                                       center + [-dx, dy]))
        self.set_verts()

    def set_verts(self):
        self.verts = self.extents.tolist() + [self.extents[0]]


class EllipticalSelection(RectangleSelection):

    def set_verts(self):
        # create an elliptical path within the extents
        max_x, max_y = np.max(self.extents, axis=0)
        min_x, min_y = np.min(self.extents, axis=0)
        a = (max_x - min_x) / 2.
        b = (max_y - min_y) / 2.
        center = [min_x + a, min_y + b]
        rad = np.arange(31) * 12 * np.pi / 180
        x = a * np.cos(rad)
        y = b * np.sin(rad)
        self.verts = np.vstack((x, y)).T + center
        self.verts = self.verts.tolist()


class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool highlights
    selected points by fading them out (i.e., reducing their alpha values).
    If your collection has alpha < 1, this tool will permanently alter them.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """
    def __init__(self, ax, collection, alpha_other=0.3, shape='Rectangle'):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, self.Npts).reshape(self.Npts, -1)

        self.selector = Selector(ax, onselect=self.onselect, shape=shape)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero([path.contains_point(xy) for xy in self.xys])[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.selector.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    #plt.ion()
    data = np.random.rand(100, 2)

    subplot_kw = dict(xlim=(0, 1), ylim=(0, 1), autoscale_on=False)
    fig, ax = plt.subplots(subplot_kw=subplot_kw)

    pts = ax.scatter(data[:, 0], data[:, 1], s=80)
    selector = SelectFromCollection(ax, pts, shape='Lasso')

    plt.show()
    raw_input('Press any key to accept selected points')
    print("Selected points:")
    print(selector.xys[selector.ind])
    selector.disconnect()

    # Block end of script so you can check that the lasso is disconnected.
    raw_input('Press any key to quit')
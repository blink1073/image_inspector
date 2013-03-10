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


class Transform(object):
    pass


#TODO: add transform tool
#TODO: allow for more than than two crossings in an intersect
#TODO: improve handling of add/sub function - better path


class Selector(AxesWidget):

    def __init__(self, ax, onselect=None, useblit=True, lineprops=None,
                 shape='Lasso'):
        AxesWidget.__init__(self, ax)

        self.useblit = useblit
        self.onselect = onselect

        if lineprops is None:
            lineprops = dict()
        lineprops.setdefault('linestyle', '--')
        self.line = Line2D([], [], **lineprops)
        self.line.set_visible(False)
        self.ax.add_line(self.line)

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

    def onkey(self, event):
        '''Update our modifiers on a key press
        '''
        if event.key in ['shift', 'alt', ' ', 'control']:
            if self.tool.active:
                self.tool.add_modifier(event.key)
            else:
                if event.key == 'shift':
                    self.mode = 'Add'
                elif event.key == 'alt':
                    self.mode = 'Subtract'
                elif event.key == 'control':
                    self.mode = 'Intersect'

    def offkey(self, event):
        '''Update our modifiers on a key release
        '''
        if event.key in self.tool.modifiers:
            self.tool.remove_modifier(event.key)

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

        if self.useblit:
            self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.line)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()

    def finalize(self):
        '''Take the appropriate action based on mode
        '''
        verts = self.tool.verts
        if self.mode in ['Add', 'Subtract', 'Intersect'] and self.ax.selection:
            try:
                verts = poly_transform(self.ax.selection, self.tool.verts,
                                        self.mode)
            except (ValueError, IndexError):
                pass
        self.onselect(verts)
        self.ax.selection = verts
        self.tool.active = False
        self.tool.verts = None
        self._patch.set_visible(False)
        self.draw_line(verts)
        self.mode = 'New'

    def set_shape(self, shape):
        if shape == 'Rectangle':
            self.tool = RectangleSelection()
        elif shape == 'Ellipse':
            self.tool = EllipticalSelection()
        elif shape == 'Lasso':
            self.tool = LassoSelection()
        elif shape == 'Polygon':
            self.tool = LassoSelection('Polygon')

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
                self.canvas.draw_idle()
        else:
            self._patch.set_visible(False)
            self.canvas.draw_idle()


class SelectorTool(object):

    def __init__(self):
        self.origin = []
        self.active = False
        self.verts = []
        self.modifiers = set()
        self.is_close = False
        self.finished = False
        self.extents = []
        self.anchor = None

    def onmove(self, event):
        pass

    def onrelease(self, event):
        pass

    def onpress(self, event):
        pass

    def finalize(self):
        pass

    def show_close(self):
        pass

    def add_modifier(self, modifier):
        self.modifiers.add(modifier)

    def remove_modifier(self, modifier):
        self.modifiers.remove(modifier)

    def start(self, event):
        self.active = True
        self.origin = (event.xdata, event.ydata)
        self.origin_pix = (event.x, event.y)
        self.verts = [self.origin]
        self.finished = False
        self.extents = np.vstack(self.verts * 4)
        self.anchor = None


class LassoSelection(SelectorTool):

    def __init__(self, type_='Lasso'):
        SelectorTool.__init__(self)
        self.type_ = type_
        self.mode = type_

    def get_mode(self):
        if ((self.type_ == 'Lasso' and 'alt' in self.modifiers) or
             (self.type_ == 'Polygon' and not 'alt' in self.modifiers)):
            return 'Polygon'
        else:
            return 'Lasso'

    def add_modifier(self, modifier):
        if modifier == 'alt' and 'alt' in self.modifiers:
            self.modifiers.remove('alt')
        else:
            self.modifiers.add(modifier)

    def remove_modifier(self, modifier):
        if not modifier == 'alt':
            self.modifiers.remove(modifier)

    def onmove(self, event):
        if not self.active or not self.origin:
            return
        if self.get_mode() == 'Polygon':
            self.verts[-1] = (event.xdata, event.ydata)
        else:
            self.verts.append((event.xdata, event.ydata))

    def onrelease(self, event):
        if not self.active:
            return
        if not self.get_mode() == 'Polygon':
            self.finalize()

    def onpress(self, event):
        if not self.verts:
            self.start(event)
            if self.get_mode() == 'Polygon':
                self.verts.append((event.xdata, event.ydata))
        elif self.get_mode() == 'Polygon':
            self.verts.append((event.xdata, event.ydata))
        if self.get_mode() == 'Polygon':
            if self.is_close and len(self.verts) > 2:
                self.finalize()

    def finalize(self):
        self.verts.append(self.origin)
        self.finished = True
        self.origin = None


class RectangleSelection(SelectorTool):

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
            if 'alt' in self.modifiers:
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

    def onrelease(self, event):
        self.finalize(event)

    def onpress(self, event):
        self.start(event)

    def finalize(self, event):
        self.finished = True

    def set_verts(self):
        self.verts = self.extents.tolist()


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


def poly_transform(old, new, transform='Add'):
    '''Add or subtract one polygon from another
    '''
    new_inner, new_outer = get_segments(old, new)
    old_inner, old_outer = get_segments(new, old)
    # find the crossover points, adjusting for path direction
    for option in ['', 'new', 'newold', 'new']:
        if 'new' in option:
            new_inner = new_inner[::-1]
            new_outer = new_outer[::-1]
        if 'old' in option:
            old_inner = old_inner[::-1]
            old_outer = old_outer[::-1]
        try:
            cross1 = segment_intersect(old_outer[-1], old_inner[0],
                                      new_outer[-1], new_inner[0])
            cross2 = segment_intersect(old_inner[-1], old_outer[0],
                                   new_inner[-1], new_outer[0])
        except ValueError:
            pass
        else:
            break
    # construct the new polygon
    if transform == 'Add':
        d1 = dist_points(new_outer[-1], cross1)
        d2 = dist_points(new_outer[-1], cross2)
        if d1 < d2:
            new_outer = new_outer[::-1]
        verts = (old_outer.tolist() + [cross1] + new_outer.tolist() + [cross2])
    elif transform == 'Subtract':
        verts = old_outer.tolist() + [cross1] + new_inner.tolist() + [cross2]
    else:
        verts = (old_inner.tolist() + [cross2] + new_inner.tolist()[::-1] +
                [cross1])
    # make it a closed polygon
    verts.append(verts[0])
    return verts


def get_segments(path1, path2):
    '''Break a polygon into the segments inside the second polygon and
outside of it
    '''
    path1 = Path(path1, closed=True)
    path2 = Path(path2, closed=True)
    locations = [path1.contains_point(xy) for xy in path2.vertices]
    crossings = np.nonzero(np.diff(locations))[0]
    if len(crossings) == 1:
        crossings = [crossings[0], -1]
    if len(crossings) != 2:
        print('Error!, incorrect number of crossings: {0}'
                .format(len(crossings)))
        raise ValueError
    else:
        verts = path2.vertices
        inner = verts[crossings[0] + 1: crossings[1] + 1]
        outer = np.vstack((verts[crossings[1] + 1:], verts[:crossings[0] + 1]))
        if locations[0]:
            inner, outer = outer, inner
        return inner, outer


def dist_points(pt1, pt2):
    '''Return the distance between two points
    '''
    x = pt1[0] - pt2[0]
    y = pt1[1] - pt2[1]
    return np.hypot(x, y)


def perpendicular(a):
    '''Return the perpendicular segment to a
    '''
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b


def segment_intersect(a1, a2, b1, b2):
    '''Find the intersection between two line segments
    '''
    da = a2 - a1
    db = b2 - b1
    dp = a1 - b1
    dap = perpendicular(da)
    denom = np.dot(dap, db)
    num = np.dot(dap, dp)
    if denom == 0:
        raise ValueError
    inter = (num / denom) * db + b1
    max_x = max(a1[0], a2[0], b1[0], b2[0])
    min_x = min(a1[0], a2[0], b1[0], b2[0])
    max_y = max(a1[1], a2[1], b1[1], b2[1])
    min_y = min(a1[1], a2[1], b1[1], b2[1])
    if (inter[0] < min_x or inter[0] > max_x or inter[1] < min_y or
            inter[1] > max_y):
        raise ValueError('Segments to not cross')
    return inter


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
    selector = SelectFromCollection(ax, pts, shape='Polygon')

    plt.show()
    raw_input('Press any key to accept selected points')
    print("Selected points:")
    print(selector.xys[selector.ind])
    selector.disconnect()

    # Block end of script so you can check that the lasso is disconnected.
    raw_input('Press any key to quit')
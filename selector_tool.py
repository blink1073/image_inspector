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
from roi import ROIToolBase


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

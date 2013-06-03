# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 21:11:25 2013

@author: silvester
"""
import numpy as np
import polygon_math

from base import CanvasToolBase
from polygon_tool import (
    RectangleSelection, LassoSelection, EllipseSelection)


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
                 maxdist=10, useblit=True, lineprops=None, shape='rectangle'):
        super(SelectionTool, self).__init__(ax, on_move=on_move,
                                on_enter=on_enter, on_release=on_release,
                                useblit=useblit)

        if on_enter is None:
            def on_enter(vertices):
                print(vertices)
        self.callback_on_enter = on_enter

        self.connect_event('key_press_event', self.onkey)
        self.connect_event('key_release_event', self.offkey)
        self.connect_custom_event('roi_changed', self.on_change)

        self.mode = 'New'
        self.verts = None
        self.modifiers = set()
        self._prev_modifiers = set()
        self._recurse = False
        self.active = False

        self.tools = [LassoSelection(ax),
                      RectangleSelection(ax, maxdist=maxdist),
                        EllipseSelection(ax, maxdist=maxdist)]

        self.shape = shape

    def onkey(self, event):
        '''Update our modifiers on a key press
        '''
        if not self.tool._busy:
            self.modifiers.add(event.key)
            self.tool._prev_line.set_visible(True)
        if event.key == 'ctrl+r':
            self.shape = 'Rectangle'
        elif event.key == 'ctrl+l':
            self.shape = 'Lasso'
        elif event.key == 'ctrl+e':
            self.shape = 'Ellipse'

    def offkey(self, event):
        '''Update our modifiers on a key release
        '''
        if not self.tool._busy and event.key in self.modifiers:
            self.modifiers.remove(event.key)

    def on_change(self, roi):
        '''Take the appropriate action based on mode
        '''
        if self._recurse:
            self._recurse = False
            return
        verts = self.tool.verts
        if self.verts is not None and self.modifiers:
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
                old = self.ax.transAxes.transform(self.verts)
                new = self.ax.transAxes.transform(self.tool.verts)
                verts = polygon_math.combine_polys(old, new, mode)
                verts = self.ax.transAxes.inverted().transform(verts).tolist()
            except (ValueError, IndexError):
                pass
            self.tool.verts = verts
            self.tool.update()
            for tool in self.tools:
                if tool.shape == 'lasso':
                    tool.verts = verts
                    self._recurse = True
                    roi.handled = True
                    tool.publish_roi()
        self.verts = verts
        self._prev_modifiers = self.modifiers
        self.modifiers = set()

    @property
    def shape(self):
        return self.tool.shape

    @shape.setter
    def shape(self, shape):
        shape = shape.lower()
        if not hasattr(self, 'tools'):
            return
        for tool in self.tools:
            if tool.shape == shape:
                tool.activate()
                self.tool = tool
                if self.verts:
                    self.tool._prev_line.set_data(zip(*self.verts))
                    self.tool._prev_line.set_visible(True)
                    self.tool.redraw()
            else:
                tool.deactivate()

    @property
    def geometry(self):
        return self.tool.verts
    def publish_roi(self):
        self.tool.publish_roi()

    def activate(self):
        self.modifiers = set()
        if hasattr(self, 'tool'):
            self.tool.activate()
        self.redraw()

    def deactivate(self):
        for tool in self.tools:
            tool.deactivate()
        self.set_visible(False)
        self.modifiers = set()
        self.verts = None
        self.redraw()

    def finalize(self):
        self.tool.finalize()
        
    @property
    def _prev_line(self):
        return self.tool._prev_line


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
    #ax.imshow(image, interpolation='nearest')
    ax.plot(np.random.random(10))
    tool = SelectionTool(ax, shape='lasso', on_release=print_geometry)
    tool.connect_event('roi_changed', roi_changed)
    plt.show()

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 20:07:11 2013

@author: silvester
"""
import numpy as np
from matplotlib.path import Path

from base import CanvasToolBase, callback_registry
from selector_tool import SelectionTool
from linetool import ThickLineTool
from point_tool import PointTool


class ROITool(CanvasToolBase):

    def __init__(self, ax, shape='rectangle',
                 on_move=None, on_enter=None, on_release=None,
                 useblit=True):
        CanvasToolBase.__init__(self, ax, on_move=on_move, on_enter=on_enter,
                            on_release=on_release, useblit=useblit)
        self.selector = SelectionTool(ax, on_release=on_release,
                                      shape='rectangle')
        self.line = ThickLineTool(ax, on_release=on_release)
        self.point = PointTool(ax, on_release=on_release)
        self.tools = [self.line, self.point, self.selector]
        self.shape = shape
        self.connect_event('roi_force', self.roi_force)
        
    def roi_force(self, roi):
        '''React to an roi setter'''
        # TODO: implement this
        import pdb; pdb.set_trace()
        pass

    def _on_key_press(self, event):
        if event.key == 'ctrl+p':
            self.activate_tool(self.point)
        elif event.key in ['ctrl+l', 'ctrl+r', 'ctrl+e']:
            self.activate_tool(self.selector)
            self.selector._on_key_press(event)
        elif event.key == 'ctrl+n':
            self.activate_tool(self.line)

    def activate_tool(self, tool):
        for item in self.tools:
            if not item is tool:
                item.deactivate()
        tool.activate()
        self._active_tool = tool
        self.redraw()

    @property
    def shape(self):
        return self.tool.shape

    @shape.setter
    def shape(self, shape):
        if shape.lower() in ['rectangle', 'ellipse', 'lasso']:
            self.activate_tool(self.selector)
            self.selector.shape = shape
        elif shape.lower() == 'line':
            self.activate_tool(self.line)
        elif shape.lower() == 'point':
            self.activate_tool(self.point)

    @property
    def geometry(self):
        return self._active_tool.geometry


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

        self.roi = ROITool(ax, on_release=self.onselect)
        self.ind = []

    def onselect(self, verts):
        verts = np.array(verts)
        self.fc[:, -1] = self.alpha_other
        if verts.size >= 5:
            path = Path(verts)
            self.ind = np.nonzero([path.contains_point(xy) for xy in self.xys])[0]
            self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.roi.disconnect()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    #plt.ion()
    data = np.random.rand(100, 2)

    def roi_changed(roi):
        print roi.shape, len(roi.geometry), type(roi.data)

    subplot_kw = dict(xlim=(0, 1), ylim=(0, 1), autoscale_on=False)
    fig, ax = plt.subplots(subplot_kw=subplot_kw)

    pts = ax.scatter(data[:, 0], data[:, 1], s=80)
    selector = SelectFromCollection(ax, pts, shape='Ellipse')
    selector.roi.connect_event('roi_changed', roi_changed)

    plt.show()
    raw_input('Press any key to accept selected points')
    print("Selected points:")
    print(selector.xys[selector.ind])
    selector.disconnect()

    # Block end of script so you can check that the lasso is disconnected.
    raw_input('Press any key to quit')

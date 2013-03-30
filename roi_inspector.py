# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:29:50 2013

@author: silvester
"""
import numpy as np

try:
    from matplotlib import lines
except ImportError:
    print("Could not import matplotlib -- skimage.viewer not available.")

from roi_tool import ROITool, CanvasToolBase


class ROIInspector(CanvasToolBase):

    def __init__(self, plot_ax, text_ax=None, use_blit=True):
        CanvasToolBase.__init__(self, plot_ax, use_blit)
        self.twin_ax = plot_ax.twinx()
        self.text_ax = text_ax
        self.connect_event('roi_changed', self.roi_changed)
        self.mode = 'histogram'
        props = dict(color='k', linewidth=1, solid_capstyle='butt')
        self._line = lines.Line2D([], [], visible=False,
                                  animated=True, **props)
        plot_ax.add_line(self._line)
        self._artists = [self._line]

    def roi_changed(self, roi):
        self.twin_ax.clear()
        if roi.shape.lower() in ['lasso', 'rectangle', 'ellipse']:
            self._line.set_visible(False)
            self.ax.clear()
            data = roi.data
            if isinstance(data, tuple):
                # TODO: use a color histogram
                pass
            else:
                nbins = min(100, np.sqrt(data.size))
                nbins = max(10, nbins)
                if self.mode == 'histogram':
                    self.ax.hist(data, bins=nbins, histtype='stepfilled')
                    self.twin_ax.hist(data, bins=nbins, color='black', normed=True, histtype='step', cumulative=True)
                else:
                    self.twin_ax.boxplot(data, 0, 'rs', 0)
            self.canvas.draw_idle()
        elif roi.shape == 'line':
            data = roi.data
            self._line.set_data((np.arange(data.size), data))
            self.set_visible(True)
            self.ax.relim()
            self.ax.autoscale_view(tight=True)
            self.redraw()
        self.display_text(roi.text)

    def display_text(self, msg):
        if not self.text_ax:
            print msg
            return
        ax = self.text_ax
        ax.clear()
        # draw the text box here
        pass
        ax.figure.canvas.draw_idle()


if __name__ == '__main__':
    np.testing.rundocs()
    import matplotlib.pyplot as plt
    from skimage import data

    image = data.camera()

    def print_geometry(geometry):
        print geometry

    def roi_changed(roi):
        print roi.shape, len(roi.geometry), type(roi.data)

    f, axes = plt.subplots(2)
    axes[0].imshow(image, interpolation='nearest')
    #ax.plot(np.random.random(10))
    tool = ROITool(axes[0], shape='line', on_release=print_geometry)
    inspector = ROIInspector(axes[1])

    plt.show()
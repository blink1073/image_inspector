# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:29:50 2013

@author: silvester
"""
import numpy as np


class ROIPlotter(object):

    def __init__(self, ax, use_blit=True):
        self.ax = ax
        self.canvas = self.ax.figure.canvas
        self.useblit = use_blit
        self.twin_ax = ax.twinx()
        self.canvas.mpl_connect('roi_changed', self.roi_changed)
        self.mode = 'histogram'
        self._line = None
        self.stat_text = ''

    def roi_changed(self, roi):
        if roi.handled:
            return
        self.twin_ax.clear()
        self.twin_ax.get_yaxis().set_visible(False)
        if roi.shape.lower() in ['lasso', 'rectangle', 'ellipse']:
            self.ax.clear()
            data = roi.data
            if roi.source_type == 'xy':
                data = data[:, 1]
            if isinstance(data, tuple):
                # TODO: use a color histogram
                pass
            elif self.mode == 'histogram':
                self.draw_histogram(data)
            elif self.mode == 'boxplot':
                self.ax.boxplot(data, 0, 'rs', 0)
            self.ax.autoscale_view(tight=True)
            self.canvas.draw_idle()
        elif roi.shape == 'line' and not roi.data is None:
            self.draw_line_profile(roi.data)
        self.stat_text = roi.stat_text

    def draw_histogram(self, data):
        nbins = min(100, np.sqrt(data.size))
        nbins = max(10, nbins)
        self.ax.hist(data, bins=nbins, histtype='stepfilled')
        self.twin_ax.hist(data, bins=nbins, color='black',
                                      normed=True, histtype='step',
                                      cumulative=True)
        self.twin_ax.get_yaxis().set_visible(True)

    def draw_line_profile(self, data):
        if not self._line or not self._line in self.ax.lines:
            self.ax.clear()
            self._line = self.ax.plot(data, 'k-')[0]
        else:
            self._line.set_xdata(np.arange(data.shape[0]))
            self._line.set_ydata(data)
        self.ax.relim()
        if self.useblit:
            self.ax.draw_artist(self._line)
        self.ax.autoscale_view(tight=True)
        self.canvas.draw_idle()


if __name__ == '__main__':
    np.testing.rundocs()
    import matplotlib.pyplot as plt
    from skimage import data
    from roi_tool import ROITool

    image = data.camera()

    def print_geometry(geometry):
        print geometry

    def roi_changed(roi):
        print roi.shape, len(roi.geometry), type(roi.data)

    f, axes = plt.subplots(2)
    axes[0].imshow(image, interpolation='nearest')
    #axes[0].plot(np.random.random(1000))
    tool = ROITool(axes[0], shape='line', on_release=print_geometry)
    inspector = ROIPlotter(axes[1])

    plt.show()
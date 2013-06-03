# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 09:29:50 2013

@author: silvester
"""
import numpy as np
from base import CanvasToolBase


class ROIPlotter(CanvasToolBase):

    def __init__(self, ax, use_blit=True):
        super(ROIPlotter, self).__init__(ax, useblit=use_blit)
        self.twin_ax = ax.twinx()
        self.ax.figure.sca(self.ax)
        self.connect_custom_event('roi_changed', self.roi_changed)
        self.mode = 'histogram'
        self._line = None
        self.stat_text = ''
        self.data = None
        self.colorbar = None
        self.busy = False

    def roi_changed(self, roi):
        if roi.handled:
            return
        if self.busy:
            return
        self.busy = True
        self.twin_ax.clear()
        self.twin_ax.get_yaxis().set_visible(False)
        if self.colorbar and not roi.shape.lower() == 'point':
            self.reset_figure()
            if roi.shape.lower == 'line':
                return
        if roi.shape.lower() in ['lasso', 'rectangle', 'ellipse']:
            self.ax.clear()
            data = roi.data
            if roi.source_type == 'xy':
                try:
                    data = data[:, 1]
                except Exception:
                    self.busy = False
                    return
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
        elif roi.shape == 'point' and not roi.data is None:
            self.ax.clear()
            try:
                im = self.ax.imshow(roi.data, interpolation='nearest', picker=True)
            except ValueError:
                self.busy = False
                return
            if not self.colorbar:
                self.colorbar = self.ax.figure.colorbar(im)
            else:
                self.colorbar.set_clim(vmin=roi.data.min(), vmax=roi.data.max())
                self.colorbar.draw_all()
            self.canvas.draw_idle()
        self.stat_text = roi.stat_text
        self.data = roi.data
        self.busy = False
    
    def reset_figure(self):
        self.ax.figure.clf()
        self.ax = self.ax.figure.add_subplot(111)
        self.twin_ax = self.ax.twinx()
        self.ax.figure.sca(self.ax)
        self.colorbar = None
        self._line = None
        self.busy = False
        self.canvas.draw_idle()
        
    def draw_histogram(self, data):
        if data is None or not data.size:
            return
        data = data[np.isfinite(data)]
        nbins = min(100, np.sqrt(data.size))
        nbins = max(10, nbins)
        try:
            self.ax.hist(data, bins=nbins, histtype='stepfilled')
        except ValueError:
            return
        vals = np.percentile(data, [2, 25, 50, 75, 98])
        for ind, val in enumerate(vals):
            if ind in [0, 4]:
                color = 'red'
            elif ind in [1, 3]:
                color = 'black'
            self.ax.axvline(x=val, ymin=0, ymax=1, color=color, linestyle='--')
        self.ax.set_xlim((vals[0], vals[-1]))
        
        #self.twin_ax.hist(data, bins=nbins, color='black',
        #                              normed=True, histtype='step',
        #                              cumulative=True)
        #self.twin_ax.get_yaxis().set_visible(True)

    def draw_line_profile(self, data):
        if not self._line or not self._line in self.ax.lines:
            self.ax.clear()
            self._line = self.ax.plot(data, 'k-')[0]
        else:
            self._line.set_xdata(np.arange(data.shape[0]))
            self._line.set_ydata(data)
        self.ax.relim()
        if self.useblit:
            self.canvas.restore_region(self.img_background)
            try:
                self.ax.draw_artist(self._line)
            except AssertionError:
                self._blit_on_draw_event()
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
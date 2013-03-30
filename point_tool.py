import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Ellipse
    LABELS_CMAP = mcolors.ListedColormap(['white', 'red', 'dodgerblue', 'gold',
                                          'greenyellow', 'blueviolet'])
except ImportError:
    print("Could not import matplotlib -- skimage.viewer not available.")

from base import CanvasToolBase
from roi import ROI


__all__ = ['PointTool']


class PointTool(CanvasToolBase):
    """Widget for painting on top of a plot.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Matplotlib axes where tool is displayed.
    overlay_shape : shape tuple
        2D shape tuple used to initialize overlay image.
    alpha : float (between [0, 1])
        Opacity of overlay
    on_move : function
        Function called whenever a control handle is moved.
        This function must accept the end points of line as the only argument.
    on_release : function
        Function called whenever the control handle is released.
    on_enter : function
        Function called whenever the "enter" key is pressed.
    rect_props : dict
        Properties for :class:`matplotlib.patches.Rectangle`. This class
        redefines defaults in :class:`matplotlib.widgets.RectangleSelector`.

    Attributes
    ----------
    overlay : array
        Overlay of painted labels displayed on top of image.
    label : int
        Current paint color.
    """
    def __init__(self, ax, radius=5, alpha=0.3, on_move=None,
                 on_release=None, on_enter=None, rect_props=None):
        super(PointTool, self).__init__(ax, on_move=on_move, on_enter=on_enter,
                                        on_release=on_release)

        props = dict(edgecolor='r', facecolor='0.7', alpha=0.5, animated=True)
        props.update(rect_props if rect_props is not None else {})

        self._point = Ellipse((0, 0), 0, 0, **props)
        self._point.set_visible(False)
        self.ax.add_patch(self._point)

        # `radius` can only be set after initializing `_point`
        self.radius = radius
        self._artists = [self._point]
        self._position = 0, 0
        self.shape = 'point'

        self.connect_event('button_press_event', self.on_mouse_press)
        self.connect_event('button_release_event', self.on_mouse_release)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, r):
        self._radius = r
        dia = 2 * r
        # translate to pixels
        p1, p2 = self.ax.transData.inverted().transform([(0, 0), (dia, dia)])
        wid = abs(p2[0] - p1[0])
        hgt = abs(p2[1] - p1[1])
        self._point.width = wid
        self._point.height = hgt

    def _on_key_press(self, event):
        if not self.active:
            return
        if event.key == 'enter':
            self.callback_on_enter(self.geometry)
            self.redraw()

    def on_mouse_press(self, event):
        if event.button != 1 or not self.ax.in_axes(event):
            return
        if not self.active:
            return
        self.update_point(event.xdata, event.ydata)

    def on_mouse_release(self, event):
        if event.button != 1 or not self.active:
            return
        self.callback_on_release(self.geometry)

    def update_point(self, x, y):
        self._point.set_visible(True)
        self._point.center = (x, y)
        self._position = (x, y)
        self.canvas.callbacks.process('roi_changed', self.roi)
        self.redraw()

    @property
    def geometry(self):
        return self._position

    @geometry.setter
    def geometry(self, pt):
        x, y = pt
        self.update_point(x, y)

    @property
    def data(self):
        if self.ax.images:
            data = self.ax.images[0].get_array()
            x, y = self._position
            return data[y, x]

    @property
    def roi(self):
        return ROI('point', self.data, self.geometry)


if __name__ == '__main__':
    np.testing.rundocs()
    import matplotlib.pyplot as plt
    from skimage import data

    image = data.camera()

    def print_geometry(geometry):
        print geometry

    f, ax = plt.subplots()
    ax.imshow(image, interpolation='nearest')
    point_tool = PointTool(ax, on_release=print_geometry)
    plt.show()

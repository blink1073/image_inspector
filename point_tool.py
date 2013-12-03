import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
except ImportError:
    print("Could not import matplotlib -- skimage.viewer not available.")

from roi import ROIToolBase


__all__ = ['PointTool']


class PointTool(ROIToolBase):
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
    def __init__(self, ax, radius=2, on_move=None,
                 on_release=None, on_enter=None, useblit=True,
                 shape_props=None):
        super(PointTool, self).__init__(ax, on_move=on_move, on_enter=on_enter,
                                        on_release=on_release, useblit=useblit)
        props = dict(edgecolor='r', facecolor='0.7', alpha=0.5, animated=True)
        props.update(shape_props if shape_props is not None else {})
        self._point = Ellipse((0, 0), 0, 0, **props)
        self._point.set_visible(False)
        self.ax.add_patch(self._point)

        # `radius` can only be set after initializing `_point`
        self._radius = radius
        self._artists = [self._point]
        self._position = 0, 0
        self.shape = 'point'

    @property
    def radii(self):
        dia = 2 * self._radius
        # translate to pixels
        # calculate asymmetry of x and y axes:
        x0, y0 = self.ax.transAxes.transform((0, 0)) # lower left in pixels
        x1, y1 = self.ax.transAxes.transform((1, 1)) # upper right in pixels
        dx = x1 - x0
        dy = abs(y1 - y0)
        maxd = max(dx, dy)
        x1, x2 = self.ax.get_xlim()
        y1, y2 = self.ax.get_ylim()
        y2 = max(y1, y2)
        width =  maxd / dx * dia / 100. * x2
        height = maxd / dy * dia / 100. * y2
        return width, height

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
        self.start(event)

    def on_mouse_release(self, event):
        if event.button != 1 or not self.active:
            return
        self.finalize()

    def update_point(self, x, y):
        self._point.width, self._point.height = self.radii
        self._point.set_visible(True)
        self._point.center = (x, y)
        self._position = (x, y)
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

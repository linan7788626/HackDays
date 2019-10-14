import numpy as np
import matplotlib.pyplot as plt


class DrawDragPoints(object):
    """
    Demonstrates a basic example of the "scaffolding" you need to efficiently
    blit drawable/draggable/deleteable artists on top of a background.
    """
    def __init__(self):
        self.fig, self.ax = self.setup_axes()
        self.xy = []
        self.tolerance = 10
        self._num_clicks = 0

        # The artist we'll be modifying...
        self.points = self.ax.scatter([], [], s=200, color='red',
                                      picker=self.tolerance, animated=True)

        connect = self.fig.canvas.mpl_connect
        connect('button_press_event', self.on_click)
        self.draw_cid = connect('draw_event', self.grab_background)

    def setup_axes(self):
        """Setup the figure/axes and plot any background artists."""
        fig, ax = plt.subplots()

        # imshow would be _much_ faster in this case, but let's deliberately
        # use something slow...
        ax.pcolormesh(np.random.random((1000, 1000)), cmap='gray')

        ax.set_title('Left click to add/drag a point\nRight-click to delete')
        return fig, ax

    def on_click(self, event):
        """Decide whether to add, delete, or drag a point."""
        # If we're using a tool on the toolbar, don't add/draw a point...
        if self.fig.canvas.toolbar._active is not None:
            return

        contains, info = self.points.contains(event)
        if contains:
            i = info['ind'][0]
            if event.button == 1:
                self.start_drag(i)
            elif event.button == 3:
                self.delete_point(i)
        else:
            self.add_point(event)

    def update(self):
        """Update the artist for any changes to self.xy."""
        self.points.set_offsets(self.xy)
        self.blit()

    def add_point(self, event):
        self.xy.append([event.xdata, event.ydata])
        self.update()

    def delete_point(self, i):
        self.xy.pop(i)
        self.update()

    def start_drag(self, i):
        """Bind mouse motion to updating a particular point."""
        self.drag_i = i
        connect = self.fig.canvas.mpl_connect
        cid1 = connect('motion_notify_event', self.drag_update)
        cid2 = connect('button_release_event', self.end_drag)
        self.drag_cids = [cid1, cid2]

    def drag_update(self, event):
        """Update a point that's being moved interactively."""
        self.xy[self.drag_i] = [event.xdata, event.ydata]
        self.update()

    def end_drag(self, event):
        """End the binding of mouse motion to a particular point."""
        for cid in self.drag_cids:
            self.fig.canvas.mpl_disconnect(cid)

    def safe_draw(self):
        """Temporarily disconnect the draw_event callback to avoid recursion"""
        canvas = self.fig.canvas
        canvas.mpl_disconnect(self.draw_cid)
        canvas.draw()
        self.draw_cid = canvas.mpl_connect('draw_event', self.grab_background)

    def grab_background(self, event=None):
        """
        When the figure is resized, hide the points, draw everything,
        and update the background.
        """
        self.points.set_visible(False)
        self.safe_draw()

        # With most backends (e.g. TkAgg), we could grab (and refresh, in
        # self.blit) self.ax.bbox instead of self.fig.bbox, but Qt4Agg, and
        # some others, requires us to update the _full_ canvas, instead.
        self.background = self.fig.canvas.copy_from_bbox(self.fig.bbox)

        self.points.set_visible(True)
        self.blit()

    def blit(self):
        """
        Efficiently update the figure, without needing to redraw the
        "background" artists.
        """
        self.fig.canvas.restore_region(self.background)
        self.ax.draw_artist(self.points)
        self.fig.canvas.blit(self.fig.bbox)

    def show(self):
        plt.show()

DrawDragPoints().show()

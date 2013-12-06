Image Inspector
===============

Brings GIMP-like selection tools to Matplotlib, plus an ROI inspector tool

Supports: Rectangle, Ellipse, Free Lasso and Polygon Lasso, Line and Point Modes


Rectangle Mode (CTRL+R)
-----------------------
- Left click and drag to drag a Rectangle
- While dragging:
   + Hold CTRL to move center to original position
   + Hold SHIFT to draw a Square

Ellipse Mode (CTRL+E)
---------------------
- Left click and drag to drag an Ellipse
- While dragging:
    + Hold CTRL to move center to original position
    + Hold SHIFT to draw a Cirlce

Lasso Mode (CTRL+L)
-------------------
- Left click and drag to free hand draw
- OR Left click and release to start a polygon segment
  Then click and release to draw a new point
- Double click at any time to close or
  Click on the start point indicator to close

Line Mode (CTRL+N)
------------------
- Left click and drag to draw a line

Point Mode (CTRL+P)
-------------------
- Left click to draw a point


Advanced Operations
-------------------
- Move a Rectangle, Ellipse, or Lasso:
        Hold SPACE and Left click within the shape to drag

- Build complex shapes
    + Hold SHIFT and create a new shape for a UNION
    + Hold CTRL and create a new shape for an INTERSECTION
    + Hold SHIFT+CTRL and create a new shape for a SUBTRACTION
    + You can switch between Rectangle, Ellipse, and Lasso at any time

# 2D-Transformation

This is the first assignment of the computer graphics course. The assignment comes with a template and my tasks were to deal with the mouse interaction and collects these points which are input by the mouse to form triangles, to draw these triangles with OpenGL function, and to generate a fractal uses these transformations based on the transformation matrices created from these triangle input.

### Implementation Workflow (My Tasks)

1. To deal with the mouse interaction and use the points which are input by the mouse to form triangles in MouseInteraction(double x, double y);
   
2. To draw triangles which have been generated by the mouse interaction in the first window in DrawTriangles();

3. To calculate the transformation matrices from the first triangle to each other triangle in AffineMatricesCalculation(Triangles in, Triangles out, Matrix matrix);

4. To use the transformation matrices calculated in CalculateMatrices() function to generate a fractal in RecursiveFractal(int k).

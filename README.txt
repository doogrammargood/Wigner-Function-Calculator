This contains python code for the Discrete Wigner Function Calculator.

The main reference is Hudson's Theorem for Finite-Dimensional Quantum Systems by David Gross, available at https://arxiv.org/abs/quant-ph/0602001

#----To Run:---#
  1.) Clone the repo with 'git clone .'
  2.) Run 'python application.py' with python version 3.5.2
  3.) Install all the packages that the program complains about. These should be
    pyqt5 version 5.14.2 (it will complain if you try to use later versions)
    pyqtgraph 0.10.0
    itertools
    pickle
    multiprocessing
  4.) After all packages are installed, 'python application.py' should run the app.


#---Current Capabilities---# (as of Monday, June 8, 2020)
  In the 'Density Matrix' tab, you can set the size of the quantum state, and choose from a list of states according to the buttons.
  Currently, the buttons '<-' '->' and the rightmost spinbox don't do anything.

  Density Matrices are considered single-particle states over a finite field.

  In the Wigner Function tab, you can see phase space in the left pane. Blue points are more positive, and yellow points are more negative.
  By hovering over a point, you'll see the identity of the point and the value of the Wigner function there in the right pane.

  If you click on a point, this will fix point1 and it will be marked with a flag. Now, hovering over a different point will display the information for this second point in point2.
  The squares which are on the line between point1 and point2 are marked with thick outlines.

  Once line is defined, you'll see a bar graph which represents the marginal probabilities of the Wigner Distribution marginalized over that line.
  The value of a line is the sum of the wigner functions of points on that line. The current line is shown in yellow in the bar graph.

  If you click on another point, the direction between the two points will be fixed. If you hover over a third point, the line will change to the parallel line through the point
  where you are hovering, and you'll see a different element of the distribution bar graph highlighted.

  Currently, Save and Edit in the toolbar don't do anything.

  Representation:
    Elements of a finite field are represented as polynomials over a prime field, modulo a Conway Polynomial (which have been obtained from the GAP system).
    When these elements are displayed, we write them as a list of integers where the ith integer in the list is the coefficient of x^i in the polynomial.
    This means that the first element of the list is the constant term of the polynomial.

    Points are represented as a pair of lists.
    Lines are represented as a triple, (a,b,c), interpreted as the line ax+by=c

    The Most Negative Point is a point where the Wigner Function is minimized
    The Total Negativity is the sum of negative Wigner Function values.

#---Planned Features---#(TODO)
  -Optimize calculating Wigner Function.
  -Allow for opening a pickled numpy matrix.
  -Implement Saving.
  -Allow hovering over elements of the bar graph to highlight the line in phase space.
  -It should say what the Conway Polynomial is somewhere.
  -Finite Field Elements should be written in Latex
  -A second view of phase space with n pxp grids (so we interpret the density matrix as describing n p-level particles rather than 1 np level particle.)
  -An editor for the density matrix:
    -Implement Clifford Covariance
    -Are there libraries to describe the clifford gates for qudits?
  -A functional tab.
    -Keep track of lines with weightings which are to be interpreted as a functional.
    -Should integrate with Clifford Covariance, so functionals can be determined by a generating line and its orbit.

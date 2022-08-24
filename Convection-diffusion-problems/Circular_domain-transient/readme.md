This script solves the transient transport problem for the circular domain already treated in the Circular_domain-CGS.ipynb file, for this reason part of the code is recycled from that script.

The script creates the matrix of coefficients using the UD interpolation for the convective method and a simple Euler implicit scheme for the time discretisation.

The system of equations is solved using the CGS method.

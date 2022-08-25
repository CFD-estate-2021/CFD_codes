The scripts in this diredctory try to solve Navier-Stokes equations for an incompressible fluid using the SIMPLE algorithm.

The domain considered is a pipe. No-slip BC are applied at the walls and either constant pressure or constant velocity BC are applied at the inlet and outlet.

The scipts were able to provide solutions for these simple settings, but the algotithm itself had important stability problems.

The case of constant velocity was very unstable due to the presence of shock waves arising from the continuity equation and the fluid incompressibility.

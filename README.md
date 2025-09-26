# Internship with Cambridge Laboratory for Scientific Computing

## Solving Euler equations in 1D & 2D
  - Using Lax Friedrichs + Richtymer fluxes and First ORder CEntered (FORCE) scheme
  - Incorporated SLIC (slope limiting) with MinBee limiter

## MHD Solvers in 1D & 2D
  - Created Godunov solvers using HLL and HLLC (2 wave and 3 wave approximations of exact Riemann problem)
  - Incorporated Van-Leer slope limiting -> MUSCL-Hancock scheme
  - Implemented Divergence Cleaning
  - Validation of sovlers using Brio-Wu (shock tube) tests in 1 and 2D
  - Further validation using Orszang Tang and Kelvin-Helmholtz tests

<table>
  <tr>
    <td>
      <img src="KelvinHelmholtz/BxGif.gif" alt="Kelvin-Helmholtz Animation" width="300"/><br/>
      <p align="center"><em>Figure 1: Animation of B_x field for Kelvin-Helmholtz Instability test using MUSCL-Hancock Scheme</em></p>
    </td>
    <td>
      <img src="OrszangTangPlots/DensityGif.gif" alt="Orszang-Tang Animation" width="300"/><br/>
      <p align="center"><em>Figure 2: Animation of Density for Orszang-Tang vortex using MUSCL-Hancock Scheme</em></p>
    </td>
  </tr>
</table>


## Linear Solver for Solov'Ev equation
  - Using C++ eigen library to solve sparse linear Solov'ev equation in 2D (axisymmetric case)

<td>
<img src="SolovEvPlots/solution.png" alt="Solov'Ev Plot" width="300"/>
<p align="left"><em>Figure 3: Sparse Linear solver solution to Solov'Ev equation</em></p>
</td> 





  

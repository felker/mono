%2D Time-Dependent RT solver
%Based on Jiang14 algorithm
%Monotonic Upwind Interpolation
%Mixed frame to O(v/c)

%Parameters
clear all;
close all;
nx = 100;
ny = 100; 
dt = 0.01; %cfl=1

nt = 160;
ntheta = 2; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lz = 1.0;
c = 1.0;
dz = lz/nx;

%Angular Discretization, Sn discrete ordinates. 
[mu, w] = angular_quad1D(ntheta);

%Spatial Discretization
%Radiation gridpoints are centered on the fluid cells
zz=linspace(0,lz,nx)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ntheta); 
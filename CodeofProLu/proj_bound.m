% Program: proj_bound.m
% To perform orthogonal projection onto set B_l,u
% see Eq. (3.31) of the notes, p. 42.
% Written by W.-S. Lu, University of Victoria.
function x = proj_bound(xi,L,U)
xw = max(xi,L);
x = min(xw,U);
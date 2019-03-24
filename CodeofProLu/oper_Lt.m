% Program: oper_Lt.m
% To perform operator LT, see Eq. (3.30) of the notes,
% p. 42.
% Written by W.-S. Lu, University of Victoria.
function [p,q] = oper_Lt(x)
[m,n] = size(x);
p = x(1:m-1,:) - x(2:m,:);
q = x(:,1:n-1) - x(:,2:n);
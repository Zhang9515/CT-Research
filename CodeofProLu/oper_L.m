% Program: oper_L.m
% To perform operator L, see Eq. (3.29) of the notes,
% p. 42.
% Written by W.-S. Lu, University of Victoria.
function x = oper_L(p,q)
[t,n] = size(p);
m = t + 1;
zn = zeros(1,n);
zm = zeros(m,1);
pe = [zn;p;zn];
qe = [zm q zm];
x = pe(2:m+1,:) + qe(:,2:n+1) - pe(1:m,:) - qe(:,1:n); 
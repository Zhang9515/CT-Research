% Program proj_pair.m
% To perform projection onto set P using the formula at the 
% bottom of p. 2425 of Reference [17].
% Written by W.-S. Lu, University of Victoria.
function [r,s] = proj_pair(p,q)
[t,n] = size(p);
m = t + 1;
p1 = p(:,1:n-1);
q1 = q(1:m-1,:);
pq1 = p1.^2 + q1.^2;
mpq = max(1,pq1);
r1 = p1./mpq;
r2 = p(:,n)./max(1,abs(p(:,n)));
r = [r1 r2];
s1 = q1./mpq;
s2 = q(m,:)./max(1,abs(q(m,:)));
s = [s1;s2];
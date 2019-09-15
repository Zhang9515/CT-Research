function C = gen_dct(n)

alpha = [sqrt(1/n) sqrt(2/n)*ones(1,n-1)];
index = (1:2:(2*n-1)) *pi / (2*n);
C = zeros(n,n);
for k = 1:n
    C(k,:) = alpha(k) * cos(index * (k-1));
    
end

end
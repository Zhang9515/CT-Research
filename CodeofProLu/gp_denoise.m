% Program: gp_denoise.m
% To remove noise from a corrupted still image by using 
% the gradient projection (GP) algorithm described in Sec. 3.2
% of the notes, see pp. 41-44.
% Reference: A. Beck and M. Teboulle, "Fast gradient-based 
% algorithms for constrained total variation image denoising 
% and deblurring problems," IEEE Trans. Image Processing, 
% vol. 18, no. 11, pp. 2419-2434, Nov. 2009.
% Input:
% xorig: original image.
% st: initial random state for noise.
% sig: standard deviation of the noise.
% lam: parameter \lambda.
% K: number of iterations to be performed.
% Output:
% x: de-noised image.
% psnr: profile of peak signal-to-noise ratio in dB.
% Written by W.-S. Lu, University of Victoria.
% Example: [x,psnr] = gp_denoise(camera256,17,0.1,0.1,300);
function [x,psnr] = gp_denoise(xorig,st,sig,lam,K)

[m,n] = size(xorig);

pk = zeros(m-1,n);       % GP median parameter
qk = zeros(m,n-1);

rng ( st ) ;
u = sig*randn(m,n);
xin = xorig + u;

k = 0;
Li = 1/(8*lam);
c = sqrt(m*n); 
psnr_before = 20*log10(c/norm(u,'fro'));
psnr = zeros(K+1,1);
psnr(1) = psnr_before;

while ( k < K )
    xpq = proj_bound(xin - lam*oper_L(pk,qk),0,1);        % Projection operation
    [pw1,qw1] = oper_Lt(xpq);       % operation of Lt
    
    pw = pk + Li*pw1;      
    qw = qk + Li*qw1;
    [pk,qk] = proj_pair(pw,qw);    % update p , q , Projection of matrix
    
    x = proj_bound(xin - lam*oper_L(pk,qk),0,1);   % update x
    k = k + 1;
    psnr(k+1) = 20*log10(c/norm(x-xorig,'fro'));
end

disp('PSNR before and after denoising:')
[psnr(1) psnr(K+1) psnr(K+1)-psnr(1)]
subplot(221)
imshow(xorig)
title('original image')
subplot(222)
imshow(xin)
title('noise contaminated image')
subplot(223)
imshow(x)
title('denoised image')
subplot(224)
plot(1:1:K+1,psnr)
grid
axis([0 K+1 psnr(1)-3 psnr(K+1)+3])
axis square
title('PSNR profile')
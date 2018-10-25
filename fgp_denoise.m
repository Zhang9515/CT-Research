% 2018/03/02 edited by ZXZ 
% Program: fgp_denoise.m
% using FISTA like tricks
% To remove noise from a corrupted still image by using the
% fast gradient projection (FGP) algorithm described in Sec. 3.2
% of the notes, see pp. 44-45.
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
% Example: [x,psnr] = fgp_denoise(camera256,17,0.1,0.1,35);
function [x,psnr] = fgp_denoise(xorig,st,sig,lam,K)

[m,n] = size(xorig);
pk = zeros(m-1,n);
rk = pk;
qk = zeros(m,n-1);
sk = qk;
randn('state',st)
u = sig*randn(m,n);
xin = xorig + u;
k = 0;
tk = 1;
Li = 1/(8*lam);
c = sqrt(m*n);
psnr_before = 20*log10(c/norm(u,'fro'));
psnr = zeros(K+1,1);
psnr(1) = psnr_before;
while  ( k < K )
    xpq = proj_bound(xin - lam*oper_L(rk,sk),0,1);
    [pw1,qw1] = oper_Lt(xpq);
    pw = rk + Li*pw1;
    qw = sk + Li*qw1;
    [pk1,qk1] = proj_pair(pw,qw);
    tk1 = 0.5*(1+sqrt(1+4*tk^2));
    ak = (tk - 1)/tk1;
    rk = pk1 + ak*(pk1-pk);
    sk = qk1 + ak*(qk1-qk);
    tk = tk1;
    pk = pk1;
    qk = qk1;
    x = proj_bound(xin - lam*oper_L(pk,qk),0,1);
    k = k + 1;
    psnr(k+1) = 20*log10(c/norm(x-xorig,'fro'));
end
disp('PSNR before and after denoising:')
[psnr(1) psnr(K+1) psnr(K+1)-psnr(1)]
% map = gray(256);
subplot(221)
imshow(xorig)
% axis square
title('original image')
subplot(222)
imshow(xin)
% axis square
title('noise contaminated image')
subplot(223)
imshow(x)
% axis square
title('denoised image')
subplot(224)
plot(1:1:K+1,psnr)
grid
axis([0 K+1 psnr(1)-3 psnr(K+1)+3])
axis square
title('PSNR profile')
xlabel('iteration')
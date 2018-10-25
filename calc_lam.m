function lam = calc_lam(I,I0,var_n,ep)
%% private function: calc_lam (by Guy Gilboa).
%% calculate scalar fidelity term of the Total-Variation
%% input: current image, original noisy image, noise variance, epsilon
%% output: lambda
%% example: lam=calc_lam(I,I0,var_n,ep)
if ~exist('ep')     % good for 256 gray-level
    ep = 1;
end
[ny,nx] = size(I); ep2 = ep^2;
% estimate derivatives
I_x = (I(:,[2:nx nx])-I(:,[1 1:nx-1]))/2;
I_y = (I([2:ny ny],:)-I([1 1:ny-1],:))/2;
I_xx = I(:,[2:nx nx])+I(:,[1 1:nx-1])-2*I;
I_yy = I([2:ny ny],:)+I([1 1:ny-1],:)-2*I;
Dp = I([2:ny ny],[2:nx nx])+I([1 1:ny-1],[1 1:nx-1]);
Dm = I([1 1:ny-1],[2:nx nx])+I([2:ny ny],[1 1:nx-1]);
I_xy = (Dp-Dm)/4;
% compute numerator and denomenator
Num = I_xx.*(ep2+I_y.^2)-2*I_x.*I_y.*I_xy+I_yy.*(ep2+I_x.^2);
Den = (ep2+I_x.^2+I_y.^2).^(3/2);
Div = Num./Den;
%%% fidelity term
lam = mean(mean(Div.*(I-I0)))./var_n;




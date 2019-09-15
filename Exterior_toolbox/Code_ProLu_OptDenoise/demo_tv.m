%%% TV demo %%
%%% By Guy Gilboa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on: [ROF92] L. Rudin, S. Osher, E. Fatemi,
% "Nonlinear Total Variation based noise removal algorithms",
% Physica D 60 259-268,1992.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simple TV denoising %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = imread('camera256.tif'); % load image
% I = double(I(11:138,11:138)); % cut a piece, convert to double
I = dicomread ( 'ActualCTImage.dcm');
winL = 0 ;    winH = 4095 ;           %  set window width
I = winL + double ( I ) / ( winH - winL ) ;        

% params
iter = 80;
% dt=0.2; eps=1;
%%% Add noise
std_n = 0.2; % Gaussian noise standard deviation
In = randn(size(I))*std_n; % White Gaussian noise
I0 = I + In; % noisy input image
% show original and noisy images

% close all
% figure(1); imshow(I,[]); title('Original')
% figure(2); imshow(I0,[]); title('Noisy image')

% denoise image by using tv for some iterations

% J = tv(I0,iter);
% figure(3); imshow(J,[]); title('Denoised image')

%%% TV denoising with (scalar) data fidelity term %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we assume the noise vaiance is known.
%% Therefore we can compute lambda (instead of fixing
%% iter). TV therfore can run automatically quite well.
%% The process is slower.
J = I0;
% params
ep_J = 1e-5; % minimum mean change in image J
lam = 0; J_old = 0;
iter = 50; dt = 5e-5; eps = 1;
var_n = std_n^2; % noise variance
i = 0;
tic
while (mean(mean(abs(J - J_old))) > ep_J) % iterate until convergence
J_old = J;
J = tv(J,iter,dt,eps,lam,I0); % scalar lam
lam = calc_lam(J,I0,var_n,eps); % update lambda (fidelity term)
disp ( mean(mean(abs(J - J_old)) ) )
end % for i
toc
figure(5); imshow(J); title('Denoised image with lambda')
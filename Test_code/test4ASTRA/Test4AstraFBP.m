% % test for ASTRA Parallel
tic
% clc
clear all
rows_and_cols = 512;
vol_geom = astra_create_vol_geom(rows_and_cols);
% 
% % parallel
det_width = 1;
det_count = rows_and_cols;
angles = linspace2(0,  pi , 180);
proj_geom = astra_create_proj_geom('parallel', det_width, det_count, angles);
% 
% % 2D Data Objects
pic = phantom(rows_and_cols);
idvol = astra_mex_data2d('create', '-vol', vol_geom, pic);
% 
idproj = astra_mex_data2d('create', '-sino', proj_geom);
% 
% % create projector
% 
projector_id = astra_create_projector('linear', proj_geom, vol_geom);

%% create forward projection
[sinogram_id, sinogram] = astra_create_sino_cuda(idvol, proj_geom, vol_geom);
% [sinogram_id, sinogram] = astra_create_sino(idvol, projector_id);
%% reconstruct
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);

% cfg = astra_struct('SART');
cfg = astra_struct('SART_CUDA');
% cfg = astra_struct('CGLS');
% cfg = astra_struct('FBP_CUDA');
% cfg = astra_struct('FBP');

% cfg.ProjectorId = projector_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;

% cfg for SART
cfg.option.ProjectionOrder = 'custom';
% cfg.option.ProjectionOrder = 'sequential';
cfg.option.ProjectionOrderList = [0:5:175 1:5:176 2:5:177 3:5:178 4:5:179];
% cfg.option.ProjectionOrderList = [0:29 30:59 60:89 90:119 120:149 150:179];
% cfg.option.ProjectionOrderList = [0:10:170 1:10:171 2:10:172 3:10:173 4:10:174 5:10:175 6:10:176 7:10:177 8:10:178 9:10:179];

% sart_id = astra_mex_algorithm('create', cfg);
% cgls_id = astra_mex_algorithm('create', cfg);
algo_id = astra_mex_algorithm('create', cfg);

astra_mex_algorithm('iterate', algo_id, 10*180);
% astra_mex_algorithm('iterate', cgls_id, 100);
% astra_mex_algorithm('run', algo_id);

V = astra_mex_data2d('get', recon_id);
figure,imshow(V, []);
toc

%% garbage disposal
astra_mex_data2d('delete', idvol, idproj, sinogram_id, recon_id);
astra_mex_projector('delete', projector_id);
astra_mex_algorithm('delete', algo_id);

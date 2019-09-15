%% Example
%% create phantom
V_exact = StandardPhantom(256);

%% create geometries and projector
proj_geom = astra_create_proj_geom('parallel', 3, 256, linspace2(0,pi,180));
vol_geom = astra_create_vol_geom(256,256);
proj_id = astra_create_projector('strip', proj_geom, vol_geom);

%% create forward projection
[sinogram_id, sinogram] = astra_create_sino(V_exact, proj_id);

%% reconstruct
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
cfg = astra_struct('FBP');
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
fbp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fbp_id);
V = astra_mex_data2d('get', recon_id);
imshow(V, [1 1.05]);

%% garbage disposal
astra_mex_data2d('delete', sinogram_id, recon_id);
astra_mex_projector('delete', proj_id);
astra_mex_algorithm('delete', fbp_id);
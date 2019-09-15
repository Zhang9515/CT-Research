% 3D FDK
clear all
tic
% load ('G:\CTcode\Code\ProjectionCone_3D\phantom512.mat')
pic = Diskphantom ( 128 ) ;

vol_geom = astra_create_vol_geom( 128 ,128, 128);
pic_id = astra_mex_data3d('create', '-vol', vol_geom, pic);


det_spacing_x = 2;
det_spacing_y = 2;
det_row_count = 256;
det_col_count = 256;
angles = linspace2(0,2*pi,360);
source_origin = 400;
origin_det = 400;

proj_geom = astra_create_proj_geom('cone',  det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles, source_origin, origin_det);

proj_id = astra_mex_data3d('create', '-proj3d', proj_geom);

[projdata_id, projdata] = astra_create_sino3d_cuda(pic, proj_geom, vol_geom);

%% reconstruct
recon_id = astra_mex_data3d('create', '-vol', vol_geom, 0);

cfg = astra_struct('FDK_CUDA');
% cfg = astra_struct('SIRT3D_CUDA');

cfg.ProjectionDataId = projdata_id;
cfg.ReconstructionDataId = recon_id;
% cfg.option.MinConstraint = 0;

algo_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', algo_id);
% astra_mex_algorithm('iterate', algo_id , 100);
% astra_mex_algorithm('get_res_norm',algo_id);

V = astra_mex_data3d('get', recon_id);
figure,imshow3Dfull (V, [ 0 1 ]);
toc

%% garbage disposal
astra_mex_data3d('delete', pic_id, proj_id);
astra_mex_algorithm('delete', algo_id);
% 
% load G:\CTcode\Code\test4ASTRA\Vsirt
% figure , plot ( 1 : 128 , squeeze( Vsirt( 64 , 64 ,  : ) )'  , 1 : 128 , squeeze( V( 64 , 64 ,  : ) )' , 1 : 128 , squeeze ( pic ( 64 , 64 ,  : ) )'  ) ;
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 128 0 2 ] ) ;




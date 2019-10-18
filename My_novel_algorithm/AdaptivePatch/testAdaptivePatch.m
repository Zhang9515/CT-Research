% 19/05/31 by ZXZ
%TV based iterative algorithm with patch prior
% set prior term as the penalty term instead of just setting as an initiate
% here patch size selected adaptively
tic
clear;
% close all;
%% CT parameter setting
% lab computer
% load 'G:\CTcode\Data\trial2D'
% server2 path
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D
% load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_HmR6
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_prior_324
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_angle72_ratio1
% save_path = '..\..\..\Data\Adaptive_patchsize_selection\PICCS_miu200_lamda0.001_p4s3sp5_CG_PreInit\' ;
% load 'E:\ZXZ\Data\trial2D_angle5'
% display parameter
displaywindow = [0 0.5] ;

% parameter define
BetaScanInt = deg2rad(5) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

pic = trial2D ; 
% pic =trial2D_prior ;
% clear trial2D
% pic = phantom(512) ;
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project

Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  

MaxP = RPic * ( 1 + 0.1 )  ;                                           
PInt = Resolution ;                      %   interval of S ( interval on the detect plain ), empircally pixel-detector ratio is related to size of image
Pdomain = - MaxP : PInt : MaxP ;                          % detective range
LP = length ( Pdomain ) ;

Ratio = 4 ;                                                           % should be smaller than 8
RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

R = zeros ( LP ,  LBeta ) ;   % create space to store fan projection
%% compute system matrix

% SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
load ..\..\..\Data\SysMatrix_fan_512
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
% clear pic
R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% R = ProjectionFan_2D ( picvector, height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% R = reshape( R , LP , LBeta ) ;
% figure,imshow(R,[])

%% iterative

Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = numel ( Display ) ;

% Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;

Display_previous = trial2D_angle72_ratio1 ;
Maxsize = 40 ;  % the actual size should be 2 * maxsize +1
size_interval = 2 ;
% patch operation parameter
% tic
[edge_map, thresh]=edge(pic, 'canny') ;
patchsize_map = bwdist(edge_map,'chessboard');
patchsize_map = (patchsize_map + double(edge_map)) ;      % add the edge
% pre-process 
patchsize_map( patchsize_map > Maxsize ) = Maxsize ;     % max control
patchsize_map = 1 + 2 * (floor((patchsize_map - 1)/size_interval+1)-1) ;   % discret by interval 2
patchsize_map = patchsize_map * 2 + 1 ;      % patch size map
patchsize_map = Img2vec_Mat2Cpp2D( patchsize_map ) ;
patch_level = Maxsize / size_interval ;
% patchsize_map = AdpativePatchSizeSelection( pic , Maxsize) ;    % reference image for patch selection can be the fbp image or the truth image
% toc
% patchsize_map_disp = Vec2img_Cpp2Mat2D( patchsize_map , height , width ) ;
disp('patch selection complete')
slidestep = [3 , 3] ;
Dicslidestep = [1 , 1] ;
sparsity = 3 ; 
% construct dictionary, because here image patches are directly used as
% atom, dictionary keeps still
HighQ_image = trial2D_prior ;
clear trial2D_prior
% construct whole range dictionary with various patch size

% versatileDictionary = cell( patch_level ,1 ) ;

% for size_index = 1 : patch_level
%     patchsize = ( ((size_index-1)*size_interval+1)*2+1 ) * [1,1] ;
%     disp(['patchsize: ', num2str(patchsize(1))])    
%     patchset_HighQ = ExtractPatch2D ( HighQ_image , patchsize , slidestep, 'NoRemoveDC' ) ;      % extract the patches from the high quality image as the atoms of the dictionary, which should be normalized later
%     versatileDictionary{size_index} = col_normalization( patchset_HighQ ) ;    % Dictionary, of which each atom has been normalized 
% end
versatileDictionary = AdaptiveDictExtract2D ( HighQ_image , patchsize_map , Dicslidestep , patch_level , size_interval ) ;
disp('dictionary complete')
%% adaptive patch operation
Display = zeros( LDisplay ,1) ;
patchset_LowQ = AdaptiveExtractPatch2D ( Display_previous , patchsize_map , slidestep, 'NoRemoveDC' ) ;    % extract the patches from the low quality image which need to be improved, store here to compute the DC later
disp('adaptive extraction complete')
%     Xintm = AdaptiveExtractPatch2D ( Display_previous, patchsize_map, 'RemoveDC' ) ;    % Xintm is set of patches which extracted from the low quality image
Xintm = patchset_LowQ ;

patchset_RemoveDC_tuple = Group_omp( versatileDictionary , Xintm, patchsize_map, [height , width],  sparsity, patch_level , size_interval) ;   % use group OMP to fit Xintm
disp('group_omp complete')
Image2D = AdaptivePatchSynthesis ( patchset_RemoveDC_tuple , patchset_LowQ, [height , width], slidestep, 'NoAddDC'  ) ;   % fuse all patches
disp('adaptive synthesis complete')
figure, imshow ( Image2D , displaywindow ) ;                     % display results
drawnow ;
disp(['PSNR:' , num2str(psnr(Image2D,pic,1))])
 toc
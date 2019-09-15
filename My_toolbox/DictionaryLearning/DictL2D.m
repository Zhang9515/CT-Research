%2019/03/21
% implement primary dictionary learning with OMP & K_SVD
clear ;
tic
%% image input
% server2 path
load 'E:\ZXZ\Data\trial2D'
load 'E:\ZXZ\Data\trial2D_prior_360'
load 'E:\ZXZ\Data\trial2D_angle5'
% lab computer
% load 'G:\CTcode\Data\trial2D'
% load 'G:\CTcode\Data\trial2D_prior'
% load 'G:\CTcode\Data\trial2D_angle5'

HighQ_image = trial2D ; 
LowQ_image = trial2D_angle5 ; 
pic = single(HighQ_image) ; 
% figure,imshow(pic,[0 0.5])
% pic = single(StandardPhantom( 256 )) ;
%% parameter define

% imgsize = size( pic ) ; 
% [ height , width ] = size ( pic ) ;              % store the size of picture
% 
% Size = [ 60 , 60 ] ;                                  % actual range
% 
% Resolution = max ( Size ) / height ;   % define the resolution of pic
% Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
% Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 
% 
% thetaint = deg2rad(5) ;                                                                 % theta unit 
% Maxtheta = deg2rad(360) ;
% thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range
% Ltheta = length ( thetaRange ) ; 
% 
% tmax = round ( Rpic * 1.1 ) ;
% t_int = 0.1 ;
% t_range = -tmax : t_int : tmax ;
% Lt = length ( t_range ) ;
% 
% % R = zeros ( Lt ,  Ltheta ) ;   
% picvector = reshape (pic, height * width , 1) ;
% R = ProjectionParallel_2D( picvector , height , width , Size , thetaRange' , t_range' ) ;     % store parallel beam projection
% R = reshape( R , Lt , Ltheta ) ;
% % figure,imshow(R,[])
% 
% Display = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
% Display = reshape( Display , height , width ) ;
% Display = Display' ;
% % figure,imshow( Display , [0 0.5])

%% patch operation
imgsize = size( pic ) ;
patchsize = [10 , 10] ; 
slidestep = [3 , 3] ;
patchset_LowQ = ExtractPatch2D ( LowQ_image , patchsize , slidestep, 'NoRemoveDC' ) ;    % extract the patches from the low quality image which need to be improved, store here to compute the DC later
patchset_HighQ = ExtractPatch2D ( HighQ_image , patchsize , slidestep, 'RemoveDC' ) ;    % extract the patches from the high quality image as the atoms of the dictionary, which should be normalized later

sparsity = 5 ; 
Dictionary = col_normalization( patchset_HighQ ) ;    % Dictionary, of which each atom has been normalized 
Xintm = ExtractPatch2D ( LowQ_image, patchsize, slidestep, 'RemoveDC' ) ;    % Xintm is set of patches which extracted from the low quality image
Alpha = omp( Dictionary , Xintm , Dictionary' * Dictionary , sparsity ) ;     % use OMP to fit Xintm
Image2D = PatchSynthesis ( Dictionary * Alpha, patchset_LowQ, patchsize, slidestep, imgsize, 'AddDC' ) ;    % fuse all patches

% figure,imshow(Image2D,[0 0.5])

rmseOri = RMSE ( LowQ_image , pic ) ;
rmseDict = RMSE ( Image2D , pic ) ;
 


















toc


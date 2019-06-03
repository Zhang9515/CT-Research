% 19/05/31 by ZXZ
%TV based iterative algorithm with patch prior
% set prior term as the penalty term instead of just setting as an initiate
% here patch size selected adaptively
tic
clear;
% close all;
%%
% lab computer
% load 'G:\CTcode\Data\trial2D'
% server2 path
load ..\Data\trial2D
load ..\Data\trial2D_prior_360
% load 'E:\ZXZ\Data\trial2D_angle5'
% display parameter
displaywindow = [0 0.5] ;

% parameter define
thetaint = deg2rad(5) ;                                                                 % theta unit 
Maxtheta = deg2rad(360) ;
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

pic = trial2D ; 
clear trial2D
% pic = phantom(512) ;
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.1 ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection
%% compute system matrix
load SysMatrix
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
clear pic
R = SysMatrix * double(picvector) ;        % generate projection with system matrix

%% iterative
Threshold = 1e-6 ;
innerTimes = 10 ;
MaxLim = 1 ; MinLim = 0 ;    % value limit
miu = 200 ;        % regularization parameter for TV
lamda1 = 0.001 * miu ;           lamda2 = 0.001 * miu ;            % relaxtion factor for split bregman( 2*miu is recommended by the classical paper)
alpha = 0.5 ;    % panelty balance between the prior and the current 
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = numel ( Display ) ;

load ..\Data\Display_previous
% Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;

Display_previous = Vec2img_Cpp2Mat2D( Display_previous , height , width ) ;
            
% patch operation parameter
Maxsize = 11 ;  % should be odd number
patchsize_map = AdpativePatchSizeSelection( Display_previous , Maxsize) ;
disp('patch selection complete')
slidestep = [3 , 3] ;
sparsity = 5 ; 
% construct dictionary, because here image patches are directly used as
% atom, dictionary keeps still
HighQ_image = trial2D_prior_360 ;
clear trial2D_prior_360
% construct whole range dictionary with various patch size
versatileDictionary = cell( (Maxsize-1)/2 ,1 ) ;
for size_index = 1 : (Maxsize-1)/2
    disp(['size: ', num2str(size_index)])
    patchsize = ( size_index * 2 + 1) * [1,1] ;
    patchset_HighQ = ExtractPatch2D ( HighQ_image , patchsize , slidestep, 'NoRemoveDC' ) ;      % extract the patches from the high quality image as the atoms of the dictionary, which should be normalized later
    versatileDictionary{size_index} = col_normalization( patchset_HighQ ) ;    % Dictionary, of which each atom has been normalized 
end
disp('dictionary complete')
%% adaptive patch operation
Display = zeros( LDisplay ,1) ;
patchset_LowQ = AdaptiveExtractPatch2D ( Display_previous , patchsize_map , 'NoRemoveDC' ) ;    % extract the patches from the low quality image which need to be improved, store here to compute the DC later
disp('adaptive extraction complete')
%     Xintm = AdaptiveExtractPatch2D ( Display_previous, patchsize_map, 'RemoveDC' ) ;    % Xintm is set of patches which extracted from the low quality image
Xintm = patchset_LowQ ;

patchset_RemoveDC_tuple = Group_omp( versatileDictionary , Xintm, patchsize_map, [height , width],  sparsity, Maxsize ) ;   % use group OMP to fit Xintm
disp('group_omp complete')
Image2D = AdaptivePatchSynthesis ( patchset_RemoveDC_tuple , patchset_LowQ, [height , width], 'NoAddDC'  ) ;   % fuse all patches
disp('adaptive synthesis complete')
figure, imshow ( Image2D , displaywindow ) ;                     % display results
drawnow ;

 toc
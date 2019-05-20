% 19/05/05 by ZXZ
%TV based iterative algorithm with patch prior
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

% SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
load SysMatrix
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
clear pic
R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% R = reshape( R , Lt , Ltheta ) ;
% figure,imshow(R,[])
% Norm = norm ( R ) ;
% Norm_pic = norm ( picvector ) ;
% % load A.mat
% % SysMatrix = A ; 
%% GPU-based projection 
% picvector = Img2vec_Mat2Cpp2D( pic ) ;
% R = ProjectionParallel_2D( picvector , height , width , Size ,thetaRange' , t_range' ) ;     % store parallel beam projection
% R2 = reshape( R2 , Lt , Ltheta ) ;
% figure,imshow(R2,[])

%% iterative
% S.Kaczmarz Method
EPS = 1e-10 ;
Threshold = 1e-6 ;
innerTimes = 10 ;
MaxLim = 1 ; MinLim = 0 ;    % value limit
miu = 200 ;        % regularization parameter for TV
lamda = 0.001 * miu ;                       % relaxtion factor for split bregman( 2*miu is recommended by the classical paper)
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = numel ( Display ) ;

% Err = R - ProjectionParallel_2D( single(Display) , height , width , Size ,thetaRange' , t_range' ) ;
% Err = R - SysMatrix * Display ;
% Residual = zeros ( Times ) ;  Residual ( 1 ) = norm ( Err ) / ( Norm ) ;      % used as stop condition
figure  % hold residual graph
load ..\Data\Display_previous
% Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;

Display_previous = Vec2img_Cpp2Mat2D( Display_previous , height , width ) ;
imshow ( Display_previous , displaywindow ) ;                     % display results
drawnow;
            
gradientMatrix_x = gradient2Dmatrix_x(height,width);
gradientMatrix_y = gradient2Dmatrix_y(height,width);

% preparation for CG
divergence_matrix = divergenceMatrix2D(height,width);
b1_CG = miu * (SysMatrix') * R ; 
iter_CG = 15 ;

% patch operation parameter
patchsize = [5 , 5] ; 
slidestep = [3 , 3] ;
sparsity = 5 ; 
% construct dictionary, because here image patches are directly used as
% atom, dictionary keeps still
HighQ_image = trial2D_prior_360 ;
clear trial2D_prior_360
patchset_HighQ = ExtractPatch2D ( HighQ_image , patchsize , slidestep, 'NoRemoveDC' ) ;    % extract the patches from the high quality image as the atoms of the dictionary, which should be normalized later
Dictionary = col_normalization( patchset_HighQ ) ;    % Dictionary, of which each atom has been normalized 

% max outloop
outeriter = 10 ;
rmse = zeros( innerTimes , outeriter ) ;                                       % judgement parameter
for outerloop = 1 : outeriter
    disp(['outerloop : ' , num2str(outerloop),'/',num2str(outeriter)])
    %% patch operation
    patchset_LowQ = ExtractPatch2D ( Display_previous , patchsize , slidestep, 'NoRemoveDC' ) ;    % extract the patches from the low quality image which need to be improved, store here to compute the DC later
       
    Xintm = ExtractPatch2D ( Display_previous, patchsize, slidestep, 'NoRemoveDC' ) ;    % Xintm is set of patches which extracted from the low quality image
    Alpha = omp( Dictionary , Xintm , Dictionary' * Dictionary , sparsity ) ;     % use OMP to fit Xintm
    Image2D = PatchSynthesis ( Dictionary * Alpha, patchset_LowQ, patchsize, slidestep, [height , width], 'NoAddDC' ) ;    % fuse all patches
    imshow ( Image2D , displaywindow ) ;                     % display results
    drawnow ;
    Display_previous = double ( Img2vec_Mat2Cpp2D( Image2D ) );
    
    % initial of inner SB iterative
    dx = gradientMatrix_x * Display_previous ; dy = gradientMatrix_y * Display_previous ; bx = zeros(LDisplay,1); by =zeros(LDisplay,1);   
    rmse(1,outerloop) = RMSE( Display_previous , picvector ) ;   % used as stop condition
    local_e = LocalError( Display_previous , Display ) ;     % initial
    IterativeTime = 1  ;      % times to iterative
    disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    rmse: ', num2str(rmse ( IterativeTime , outerloop))]) ;
    %% split bregman to solve tv-based problem
    while ( IterativeTime <= innerTimes && local_e > Threshold )            % end condition of loop

                 b_CG = b1_CG + lamda * ( gradientMatrix_x * ( dx - bx ) + gradientMatrix_y * ( dy - by ) ) ;

                 Display = cgls4TV ( SysMatrix, divergence_matrix , b_CG , iter_CG , miu , lamda, Display_previous) ;

                 Display ( Display < MinLim ) = MinLim ;       Display ( Display > MaxLim ) = MaxLim ;   % non-negation constraint

                 dx = soft_threshold( gradientMatrix_x * Display + bx , 1/lamda);         % split bregman update
                 dy = soft_threshold( gradientMatrix_y * Display + by , 1/lamda);
                 bx = bx + dx - gradientMatrix_x * Display ; 
                 by = by + dy - gradientMatrix_y * Display ; 

                 IterativeTime = IterativeTime + 1 ;
                 rmse ( IterativeTime ,outerloop) = RMSE ( Display , picvector) ;      % compute error
                 local_e = LocalError( Display , Display_previous ) ;
                 loss = norm(gradientMatrix_x * Display,1) + norm(gradientMatrix_y * Display,1) + miu * norm(SysMatrix * Display - R ,2) / 2 ;     % objective function
                 disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    RMSE: ', num2str(rmse ( IterativeTime ,outerloop)),';   |    local_e: ', num2str(local_e), ';   |    Loss: ', num2str(loss)]) ;
    %     plot ( 2 : IterativeTime , rmse ( 2  : IterativeTime ) ) ;
    %     ylim ( [ 0 , ( 10 * rmse ( IterativeTime ) ) ] ) ;
    %     drawnow ;           

                Display_previous = Display ;
    end
    save ..
    Display_previous = Vec2img_Cpp2Mat2D( Display_previous , height , width ) ;
    imshow ( Display_previous , displaywindow ) ;                     % display results
    drawnow;
    clear dx dy bx by;
end
% figure , imshow ( Display_previous , displaywindow ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;


 toc
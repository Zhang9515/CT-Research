% 19/10/09 by ZXZ
%TV based iterative algorithm
% fan beam
tic
clear;
% close all;
%% CT parameter setting
% lab computer
% load 'G:\CTcode\Data\trial2D'
% server2 path
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_prior_324
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_angle72_ratio1
save_path = '..\..\..\Data\Adaptive_patchsize_selection\PICCS_miu200_lamda0.001_p4s3sp5_CG_PreInit\' ;
% load 'E:\ZXZ\Data\trial2D_angle5'
% display parameter
displaywindow = [0 0.5] ;

% parameter define
BetaScanInt = deg2rad(5) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

pic = trial2D ; 
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
clear pic
R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% R = ProjectionFan_2D ( picvector, height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% R = reshape( R , LP , LBeta ) ;
% figure,imshow(R,[])
% Norm = norm ( R ) ;
% Norm_pic = norm ( picvector ) ;
% % load A.mat
% % SysMatrix = A ; 

%% iterative
% S.Kaczmarz Method
Threshold = 1e-6 ;
MaxLim = 1 ; MinLim = 0 ;    % value limit
miu = 200 ;        % regularization parameter for TV
lamda1 = 0.001 * miu ;           lamda2 = 0.001 * miu ;            % relaxtion factor for split bregman( 2*miu is recommended by the classical paper)
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = numel ( Display ) ;

% Err = R - ProjectionParallel_2D( single(Display) , height , width , Size ,thetaRange' , t_range' ) ;
% Err = R - SysMatrix * Display ;
% Residual = zeros ( Times ) ;  Residual ( 1 ) = norm ( Err ) / ( Norm ) ;      % used as stop condition
figure  % hold residual graph
% Display_previous = trial2D_angle72_ratio1 ;
Display_previous = zeros(LDisplay,1) ;
% Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;

drawnow;
            
gradientMatrix_x = gradient2Dmatrix_x(height,width);
gradientMatrix_y = gradient2Dmatrix_y(height,width);

% preparation for CG
divergence_matrix = divergenceMatrix2D(height,width);
b1_CG = miu * (SysMatrix') * R ; 
% b1_CG = miu * BackprojectionFan2D( single(R) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;
iter_CG = 80 ;
% img_init = zeros(LDisplay,1) ; 

% set loop times
outeriter = 20 ;
innerTimes = 10 ;

rmse = zeros( innerTimes , outeriter ) ;                                       % judgement parameter
PSNR = zeros( innerTimes , outeriter ) ;
% texture part: 280:360,120:200
% plain part : 140:220,320:400

for outerloop = 1 : outeriter
    disp(['outerloop : ' , num2str(outerloop),'/',num2str(outeriter)])
    %% patch operation
    Display = double ( Img2vec_Mat2Cpp2D( Display_previous ) ) ;      % set the result of last iterate as the initiate value of current iterate
    
    % initial of inner SB iterative, all parameter should be initialized
    % based on zero input
    dx1 = zeros(LDisplay,1); dy1 = zeros(LDisplay,1); bx1 = zeros(LDisplay,1); by1 =zeros(LDisplay,1);   
    
    % initial of range limitation parameter

    local_e = 100 ;     % initial
    IterativeTime = 1  ;      % times to iterative
    
%     Display = zeros(LDisplay,1) ;
    %% split bregman to solve tv-based problem
    while ( IterativeTime <= innerTimes && local_e > Threshold )            % end condition of loop
                 Display_previous = Display ;
                 % introduce the range limitation into the iterative
                 % framework

                 Display_det1 = 0 ; Display_det2 = 0 ;
                 
                 b_CG = b1_CG  ;  %+ lamda1 * ( gradientMatrix_x * ( dx1 - bx1 ) + gradientMatrix_y * ( dy1 - by1 ) )
                  
                 % here are two choices: 1. using the previous result; 2. using zero initialization
                 % To solve Ax = b , the paramter matrix A is already a
                 % semidefinite full-rank matrix, so there is no need to
                 % use norm equation
                 systemfunc = @(vec,tt) miu * (SysMatrix') * ( SysMatrix * vec )  ;   %+ lamda1 * divergence_matrix * vec
%                  systemfunc = @(vec,tt) miu * BackprojectionFan2D( single ( ProjectionFan_2D ( single(vec), height, width, Size, BetaScanRange', Pdomain', RScan ) ) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ...
%                      + (lamda1 + lamda2) * divergence_matrix * vec ;
%                  [Display,flag,resNE,iter] = cgls(systemfunc, b_CG, 0, 1e-6, iter_CG, true, Display_previous) ;
                 [Display,flag,resNE,iter] = cg(systemfunc, b_CG, 1e-6, iter_CG, true, Display_previous, picvector) ;
                 
%                  Display = cg4TV ( SysMatrix, divergence_matrix , Display_det1, Display_det2, b_CG , iter_CG , miu , lamda1 + lamda2 , Display_previous ) ;        
%                  Display = cgls4TV ( SysMatrix, divergence_matrix , b_CG , iter_CG , miu , lamda1 + lamda2 , Display_previous) ;        
%                  miu1 = max(0 , miu1 - xig1 * (Display - MinLim)) ; miu2 = max( 0 , miu2 - xig2 * (-Display + MaxLim ) ) ;
%                  xig1 = xig1 * rate1 ; xig2 = xig2 * rate2 ;
                 Display ( Display < MinLim ) = MinLim ;       Display ( Display > MaxLim ) = MaxLim ;   % non-negation constraint
                 
                 dx1 = soft_threshold( gradientMatrix_x * Display + bx1 , 1/(lamda1+eps));         % split bregman update
                 dy1 = soft_threshold( gradientMatrix_y * Display + by1 , 1/(lamda1+eps));
                 bx1 = Cut( gradientMatrix_x * Display + bx1 , 1/(lamda1+eps) ) ; 
                 by1 = Cut( gradientMatrix_y * Display + by1 , 1/(lamda1+eps) ) ; 
             
                 rmse ( IterativeTime ,outerloop) = RMSE ( Display , picvector) ;      % compute error
                 PSNR( IterativeTime ,outerloop) = psnr ( Display , double(picvector) , 1) ; 
                 local_e = LocalError( Display , Display_previous ) ;
                 % objective function ( which is different from the previous)
                 loss = (norm(gradientMatrix_x * Display,1) + norm(gradientMatrix_y * Display,1)) + ...
                 + 0.5 * miu * norm(ProjectionFan_2D ( single(Display), height, width, Size, BetaScanRange', Pdomain', RScan ) - R ,2) ;     
                 disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    RMSE: ', num2str(rmse ( IterativeTime ,outerloop)), ';   |    psnr: ', num2str(PSNR ( IterativeTime , outerloop)), ';   |    local_e: ', num2str(local_e), ';   |    Loss: ', num2str(loss)]) ;
                 disp( ['SplitBregman_Constraint_X1: ',  num2str(norm(dx1-gradientMatrix_x * Display,2)), '  Constraint_Y1: ', num2str(norm(dy1-gradientMatrix_y * Display,2))] )
    %     plot ( 2 : IterativeTime , rmse ( 2  : IterativeTime ) ) ;
    %     ylim ( [ 0 , ( 10 * rmse ( IterativeTime ) ) ] ) ;
    %     drawnow ;           

                IterativeTime = IterativeTime + 1 ;
                test_display = Vec2img_Cpp2Mat2D( Display , height , width ) ;
                imshow ( test_display , displaywindow ) ;                     % display results
                drawnow;
    end

    Display = Vec2img_Cpp2Mat2D( Display , height , width ) ;
    imshow ( Display , displaywindow ) ;                     % display results
    drawnow;
%     save_path_pic = strcat(save_path,num2str(outerloop)) ;
%     save( save_path_pic, 'Display') ;
    Display_previous = Display ;
end
% save_path_rmse = strcat(save_path , 'rmse') ;
% save( save_path_rmse , 'rmse' );
% save_path_psnr = strcat(save_path , 'PSNR') ;
% save( save_path_psnr , 'PSNR' );
figure , imshow ( Display , displaywindow ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;
%% save image code
% figure,imshow ( test_display, displaywindow,'border', 'tight','initialmagnification','fit') ;
% set (gcf,'Position',[0,0,height,width]);   
% print -djpeg -r600 ..\..\..\Data\Adaptive_patchsize_selection\global_gold_iter5_SB_out1   % grey image -r600, color image -r300
% figure,imshow ( p30_iter2_zoom_plain_out3, displaywindow,'border', 'tight','initialmagnification','fit') ;
% set (gcf,'Position',[0,0,height,width]);   
% print -djpeg -r600 ..\..\..\Data\Adaptive_patchsize_selection\global_gold_iter5_zoom_plain_SB_out1    % grey image -r600, color image -r300
% figure,imshow ( p30_iter2_zoom_texture_out3, displaywindow,'border', 'tight','initialmagnification','fit') ;
% set (gcf,'Position',[0,0,height,width]);   
% print -djpeg -r600 ..\..\..\Data\Adaptive_patchsize_selection\global_gold_iter5_zoom_texture_SB_out1    % grey image -r600, color image -r300
%% process code
% trial2D_zoom_plain = trial2D(280 :360 ,120 :200) ;
% trial2D_zoom_texture = trial2D(140 :220 ,320 :400 ) ;
% p30_iter2_zoom_plain_out3 = test_display(280 :360 ,120 :200 ) ;
% p30_iter2_zoom_texture_out3 = test_display(140 :220 ,320 :400) ;
% psnr( test_display,trial2D)
% psnr(p30_iter2_zoom_plain_out3,trial2D_zoom_plain)
% psnr(p30_iter2_zoom_texture_out3,trial2D_zoom_texture)
% ssim( test_display,trial2D) 
% ssim(p30_iter2_zoom_plain_out3,trial2D_zoom_plain)
% ssim(p30_iter2_zoom_texture_out3,trial2D_zoom_texture)
% rRMSE( test_display,trial2D)
% rRMSE(p30_iter2_zoom_plain_out3,trial2D_zoom_plain)
% rRMSE(p30_iter2_zoom_texture_out3,trial2D_zoom_texture)


 toc
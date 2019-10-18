% 19/10/17 by ZXZ
%TV based iterative algorithm with patch prior
% set prior term as the penalty term instead of just setting as an initiate
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

%% iterative
% S.Kaczmarz Method
soft_threshold_im = 0.001 ;
soft_threshold_pr = 0.01 ;
Threshold = 1e-6 ;
MaxLim = 1 ; MinLim = 0 ;    % value limit
% tradition parameter
% miu = 50 ;        % regularization parameter for TV
% alpha = 0.5 ;    % panelty balance between the prior and the current 
% lamda1 = ( 0.02 * alpha ) * miu  ;           lamda2 = ( 0.05 * ( 1 - alpha ) ) * miu ;            % relaxtion factor for split bregman( 2*miu is recommended by the classical paper)
% novel parameter
alpha = 0.5 ;
lamda1 = alpha / soft_threshold_im ;           lamda2 = ( 1 - alpha ) / soft_threshold_pr ;
miu = 100 * lamda1 ;
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = numel ( Display ) ;
% set loop times
innerTimes = 50 ;
iter_CG = 30 ;

% Err = R - ProjectionParallel_2D( single(Display) , height , width , Size ,thetaRange' , t_range' ) ;
% Err = R - SysMatrix * Display ;
% Residual = zeros ( Times ) ;  Residual ( 1 ) = norm ( Err ) / ( Norm ) ;      % used as stop condition
figure  % hold residual graph
Display_previous = trial2D_angle72_ratio1 ;
% Display_previous = zeros(LDisplay,1) ;
% Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;

drawnow;
            
gradientMatrix_x = gradient2Dmatrix_x(height,width);
gradientMatrix_y = gradient2Dmatrix_y(height,width);

% preparation for CG
divergence_matrix = divergenceMatrix2D(height,width);
b1_CG = miu * (SysMatrix') * R ; 
% b1_CG = miu * BackprojectionFan2D( single(R) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;

rmse = zeros( innerTimes ,1 ) ;                                       % judgement parameter
PSNR = zeros( innerTimes ,1) ;
% texture part: 280:360,120:200
% plain part : 140:220,320:400

% Display_previous = zeros( 512 , 512) ;

%% patch operation
Display = double ( Img2vec_Mat2Cpp2D( Display_previous ) ) ;      % set the result of last iterate as the initiate value of current iterate

Display_prior = double ( Img2vec_Mat2Cpp2D( trial2D_prior ) ) ;
% Display_prior = double(trial2D_prior) ;

% initial of inner SB iterative, all parameter should be initialized
% based on zero input
%     dx1 = zeros(LDisplay,1); dy1 = zeros(LDisplay,1); bx1 = zeros(LDisplay,1); by1 =zeros(LDisplay,1);   
dx1 = gradientMatrix_x * Display ; dy1 = gradientMatrix_y * Display ; bx1 = zeros(LDisplay,1); by1 =zeros(LDisplay,1);   
%     dx2 = gradientMatrix_x * (-Display_prior) ; dy2 = gradientMatrix_y * (-Display_prior) ; bx2 = zeros(LDisplay,1); by2 =zeros(LDisplay,1);   
dx2 = gradientMatrix_x * (Display - Display_prior) ; dy2 = gradientMatrix_y * (Display-Display_prior) ; bx2 = zeros(LDisplay,1); by2 =zeros(LDisplay,1);  
%     dx2 = zeros(LDisplay,1) ; dy2 = zeros(LDisplay,1) ; bx2 = zeros(LDisplay,1); by2 =zeros(LDisplay,1); 
% initial of range limitation parameter
% parameter initiation in the range limitation ( min < x < max)
miu1 = zeros(LDisplay,1) ;
miu2 = zeros(LDisplay,1) ;
xig1 = 0.1 * 2^0;
xig2 = 0.1 * 2^0;
rate1 = 2 ; 
rate2 = 2 ;

local_e = 100 ;     % initial
IterativeTime = 1  ;      % times to iterative

%     Display = zeros(LDisplay,1) ;
%% split bregman to solve tv-based problem
while ( IterativeTime <= innerTimes && local_e > Threshold )            % end condition of loop
             Display_previous = Display ;
             % introduce the range limitation into the iterative
             % framework
%                  Display_det1 = max( 0 , miu1 - xig1 * ( Display_previous - MinLim ) ) ;  Display_det1 = Display_det1 ./ ( Display_det1 + eps ) ;
%                  Display_det2 = max( 0 , miu2 - xig2 * ( -Display_previous + MaxLim ) ) ;  Display_det2 = Display_det2 ./ ( Display_det2 + eps ) ;
%                  Display_det1 = xig1 * sparse( 1:LDisplay, 1:LDisplay, Display_det1, LDisplay , LDisplay ) ;
%                  Display_det2 = xig2 * sparse( 1:LDisplay, 1:LDisplay, Display_det2, LDisplay , LDisplay ) ;
             Display_det1 = 0 ; Display_det2 = 0 ;

             b_CG = b1_CG + lamda1 * ( gradientMatrix_x' * ( dx1 - bx1 ) + gradientMatrix_y' * ( dy1 - by1 ) ) ... 
             + lamda2 * ( gradientMatrix_x' * ( dx2 + gradientMatrix_x * Display_prior - bx2 ) + gradientMatrix_y' * ( dy2 + gradientMatrix_y * Display_prior - by2 ) ) ...
             + Display_det1 * ( xig1 * MinLim + miu1 ) + Display_det2 * ( xig2 * MaxLim - miu2 ) ;

             % here are two choices: 1. using the previous result; 2. using zero initialization
             % To solve Ax = b , the paramter matrix A is already a
             % semidefinite full-rank matrix, so there is no need to
             % use norm equation
             systemfunc = @(vec,tt) miu * (SysMatrix') * ( SysMatrix * vec ) + (lamda1 + lamda2) * divergence_matrix * vec ;
%                  systemfunc = @(vec,tt) miu * BackprojectionFan2D( single ( ProjectionFan_2D ( single(vec), height, width, Size, BetaScanRange', Pdomain', RScan ) ) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ...
%                      + (lamda1 + lamda2) * divergence_matrix * vec ;
%                  [Display,flag,resNE,iter] = cgls(systemfunc, b_CG, 0, 1e-6, iter_CG, true, Display_previous) ;
             [Display,flag,resNE,iter] = cg(systemfunc, b_CG, 1e-6, iter_CG, true, Display_previous, picvector) ;

%                  Display = cg4TV ( SysMatrix, divergence_matrix , Display_det1, Display_det2, b_CG , iter_CG , miu , lamda1 + lamda2 , Display_previous ) ;        
%                  Display = cgls4TV ( SysMatrix, divergence_matrix , b_CG , iter_CG , miu , lamda1 + lamda2 , Display_previous) ;        
%                  miu1 = max(0 , miu1 - xig1 * (Display - MinLim)) ; miu2 = max( 0 , miu2 - xig2 * (-Display + MaxLim ) ) ;
%                  xig1 = xig1 * rate1 ; xig2 = xig2 * rate2 ;
             Display ( Display < MinLim ) = MinLim ;       Display ( Display > MaxLim ) = MaxLim ;   % non-negation constraint

             Substract_Display = Display - Display_prior ;
             dx1 = soft_threshold( gradientMatrix_x * Display + bx1 , alpha/lamda1);         % split bregman update
             dy1 = soft_threshold( gradientMatrix_y * Display + by1 , alpha/lamda1);
             dx2 = soft_threshold( gradientMatrix_x * Substract_Display + bx2 , (1-alpha)/lamda2);         % split bregman update
             dy2 = soft_threshold( gradientMatrix_y * Substract_Display + by2 , (1-alpha)/lamda2);
             bx1 = Cut( gradientMatrix_x * Display + bx1 , alpha/lamda1 ) ; 
             by1 = Cut( gradientMatrix_y * Display + by1 , alpha/lamda1 ) ; 
             bx2 = Cut( gradientMatrix_x * Display + bx2 , (1-alpha)/lamda2 ) ; 
             by2 = Cut( gradientMatrix_y * Display + by2 , (1-alpha)/lamda2 ) ; 

             rmse ( IterativeTime) = RMSE ( Display , picvector) ;      % compute error
             PSNR( IterativeTime) = psnr ( Display , double(picvector) , 1) ; 
             local_e = LocalError( Display , Display_previous ) ;
             % objective function ( which is different from the previous)
%                  loss = alpha * (norm(gradientMatrix_x * Display,1) + norm(gradientMatrix_y * Display,1)) + ( 1 - alpha ) * (norm(gradientMatrix_x * Substract_Display,1) + norm(gradientMatrix_y * Substract_Display,1))...
%                  + 0.5 * miu * norm(ProjectionFan_2D ( single(Display), height, width, Size, BetaScanRange', Pdomain', RScan ) - R ,2) ;     
             loss = alpha * (norm(gradientMatrix_x * Display,1) + norm(gradientMatrix_y * Display,1)) + ( 1 - alpha ) * (norm(gradientMatrix_x * Substract_Display,1) + norm(gradientMatrix_y * Substract_Display,1))...
             + 0.5 * miu * norm(SysMatrix * double(Display) - R ,2) ;    
             disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    RMSE: ', num2str(rmse ( IterativeTime)), ';   |    psnr: ', num2str(PSNR ( IterativeTime)), ';   |    local_e: ', num2str(local_e), ';   |    Loss: ', num2str(loss)]) ;
             disp( ['SplitBregman_Constraint_X1: ',  num2str(norm(dx1-gradientMatrix_x * Display,2)), '  Constraint_Y1: ', num2str(norm(dy1-gradientMatrix_y * Display,2)), ...
                 '  Constraint_X2: ', num2str(norm(dx2-gradientMatrix_x * Substract_Display,2)), '  Constraint_Y2: ' , num2str(norm(dy2-gradientMatrix_y * Substract_Display,2))] )
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

save_path_rmse = strcat(save_path , 'rmse') ;
save( save_path_rmse , 'rmse' );
save_path_psnr = strcat(save_path , 'PSNR') ;
save( save_path_psnr , 'PSNR' );
% figure , imshow ( Display_previous , displaywindow ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;
%% save image code
figure,imshow ( test_display, displaywindow,'border', 'tight','initialmagnification','fit') ;
set (gcf,'Position',[0,0,height,width]);   
print -djpeg -r600 ..\..\..\Data\Adaptive_patchsize_selection\test_display   % grey image -r600, color image -r300
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
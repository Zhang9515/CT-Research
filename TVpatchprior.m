% 19/05/05 by ZXZ
%TV based iterative algorithm with patch prior
tic
clear;
% close all;
%%
% lab computer
% load 'G:\CTcode\Data\trial2D'
% server2 path
load 'E:\ZXZ\Data\trial2D'

% parameter define
thetaint = deg2rad(5) ;                                                                 % theta unit 
Maxtheta = deg2rad(360) ;
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

% pic = trial2D ; 
pic = phantom(256) ;
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
% R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% R = reshape( R , Lt , Ltheta ) ;
% figure,imshow(R,[])
% Norm = norm ( R ) ;
Norm_pic = norm ( picvector ) ;
% % load A.mat
% % SysMatrix = A ; 
%% GPU-based projection 
% picvector = Img2vec_Mat2Cpp2D( pic ) ;
R = ProjectionParallel_2D( picvector , height , width , Size ,thetaRange' , t_range' ) ;     % store parallel beam projection
% R2 = reshape( R2 , Lt , Ltheta ) ;
% figure,imshow(R2,[])

%% iterative
% S.Kaczmarz Method
EPS = 1e-10 ;
Threshold = 1e-6 ;
Times = 100 ;
IterativeTime = 1  ;      % times to iterative
miu = 0.01 ;        % regularization parameter for TV
lamda = 2 * miu ;                       % relaxtion factor for split bregman( 2*miu is recommended by the classical paper)
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
LDisplay = size ( Display ) ;
ME = zeros( 1 ,  Times ) ;                                       % judgement parameter

% Err = R - ProjectionParallel_2D( single(Display) , height , width , Size ,thetaRange' , t_range' ) ;
% Err = R - SysMatrix * Display ;
% Residual = zeros ( Times ) ;  Residual ( 1 ) = norm ( Err ) / ( Norm ) ;      % used as stop condition
figure  % hold residual graph

dx = zeros(LDisplay); dy = zeros(LDisplay); bx = zeros(LDisplay); by =zeros(LDisplay);
Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
ME(1) = norm ( Display_previous - picvector ) /  ( Norm_pic ) ;   % used as stop condition

gradientMatrix_x = gradient2Dmatrix_x(height,width);
gradientMatrix_y = gradient2Dmatrix_y(height,width);
% preparation for CG
divergence_matrix = divergenceMatrix2D(height,width);
b1_CG = miu * (SysMatrix') * R ; 
iter_CG = 5 ;

while ( IterativeTime <= Times && ME ( IterativeTime ) > Threshold )            % end condition of loop
                      
             b_CG = b1_CG + lamda * ( gradientMatrix_x * ( dx - bx ) + gradientMatrix_y * ( dy - by ) ) ;
             
             Display = cgls4TV ( SysMatrix, divergence_matrix , b_CG , iter_CG , miu , lamda, Display_previous) ;
            
             Display ( Display < 0 ) = 0 ;       % non-negation constraint
             
             dx = soft_threshold( gradientMatrix_x * Display + bx , 1/lamda);         % split bregman update
             dy = soft_threshold( gradientMatrix_y * Display + by , 1/lamda);
             bx = bx + dx - gradientMatrix_x * Display ; 
             by = by + dy - gradientMatrix_y * Display ; 
             
             IterativeTime = IterativeTime + 1 ;
             ME ( IterativeTime ) = norm ( Display - picvector ) /  ( Norm_pic ) ;      % compute error
             local_e = norm(Display - Display_previous,2) ;
             disp ( ['IterativeTime: ', num2str(IterativeTime), '; ME: ', num2str(ME ( IterativeTime )),'; local_e: ', num2str(local_e)]) ;
%     plot ( 2 : IterativeTime , ME ( 2  : IterativeTime ) ) ;
%     ylim ( [ 0 , ( 10 * ME ( IterativeTime ) ) ] ) ;
%     drawnow ;           
            Display_med = Vec2img_Cpp2Mat2D( Display , height , width ) ;
            imshow ( Display_med , [ 0  0.5 ] ) ;                     % display results
            drawnow;
            Display_previous = Display ;
end

Display = Vec2img_Cpp2Mat2D( Display , height , width ) ;

figure , imshow ( Display , [ 0  0.5 ] ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;


 toc
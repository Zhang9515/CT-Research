%% SIRT 17/12/16 ZXZ
% update at 19/05/04
% actually what i do is SIRT , which is different from SART in the term of projection line. 
tic
clear all;
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

% pic = phantom ( 257 ) ;
% pic = dicomread ( 'ActualCTImage.dcm');
% winL = 0 ;    winH = 4095 ;           %  set window width
% pic = winL + double ( pic ) / ( winH - winL ) ;    

pic = trial2D ; 
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
%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix

% SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
load SysMatrix
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
% R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% R = reshape( R , Lt , Ltheta ) ;
% figure,imshow(R,[])
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
Lamda = 1 ;                       % relaxtion factor  SIRT
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
MinLim = 0 ; MaxLim = 1;     % range constraint of solution
% Display = ones ( 1  , height * width) ;          % store the reconstruction for MART
LDisplay = size ( Display ) ;
rmse = zeros( 1 ,  Times ) ;                                       % judgement parameter
Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
Err = R - SysMatrix * Display_previous ;
local_e = LocalError( Display , Display_previous ) ;
figure  % hold residual graph

% compute the row/column sum of system matrix
Rowsum = sum( SysMatrix , 2 ) ;
Colsum = sum ( SysMatrix ) ;

while ( IterativeTime <= Times && local_e > Threshold )            % end condition of loop
%              
           
             Display = Display_previous + Lamda * SysMatrix' * ( Err ./ ( Rowsum + EPS ) ) ./ ( Colsum + EPS )' ; 
             
             Display ( Display < MinLim ) = MinLim ;       Display ( Display > MaxLim ) = MaxLim ;   % non-negation constraint
%     Dismean = mean ( Display ) ;
            rmse ( IterativeTime ) = RMSE( Display , picvector ) ;      % compute error
            local_e = LocalError( Display , Display_previous ) ;
            disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    RMSE: ', num2str( rmse(IterativeTime )),';   |    local_e: ', num2str(local_e) ]) ;
            IterativeTime = IterativeTime + 1 ;
             
            Display_med = Vec2img_Cpp2Mat2D( Display , height , width ) ;
            imshow ( Display_med , [ 0  0.5 ] ) ;                     % display results
            drawnow;
            Display_previous = Display ;
            Err = R - SysMatrix * Display ;
%     Rmed = ProjectionParallel_2D( single(Display) , height , width , Size , deg2rad(thetaRange') , t_range' ) ;
%     Residual ( IterativeTime ) = norm ( R - Rmed ) / ( Norm * height * width ) ;  
    
%     plot ( 2 : IterativeTime , Residual ( 2  : IterativeTime ) ) ;
%     ylim ( [ 0 , ( 10 * Residual ( IterativeTime ) ) ] ) ;
%     drawnow ; 
    
end

Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;

figure , imshow ( Display , [ 0  0.5 ] ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;


 toc
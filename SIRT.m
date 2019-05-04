%% SIRT 17/12/16 ZXZ
% actually what i do is SIRT , which is different from SART in the term of projection line. 
tic
clear all;
% close all;
%%
% lab computer
load 'G:\CTcode\Data\trial2D'

% parameter define
thetaint = 5 ;                                                                 % theta unit 
Maxtheta = 360 ;
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

% pic = phantom ( 257 ) ;
% pic = dicomread ( 'ActualCTImage.dcm');
% winL = 0 ;    winH = 4095 ;           %  set window width
% pic = winL + double ( pic ) / ( winH - winL ) ;        
pic = trial2D ; 

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

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
picvector = reshape ( flipud( pic )' , height * width , 1  ) ;  % original image
R = SysMatrix * picvector ;        % generate projection with system matrix
Norm = norm ( R ) ;
Norm_pic = norm ( picvector ) ;
% load A.mat
% SysMatrix = A ; 
%% iterative
% S.Kaczmarz Method
EPS = 1e-10 ;
Times = 1000 ;
IterativeTime = 1  ;      % times to iterative
Lamda = 1 ;                       % relaxtion factor  SIRT
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
% Display = ones ( 1  , height * width) ;          % store the reconstruction for MART
LDisplay = size ( Display ) ;
ME = zeros( 1 ,  Times ) ;                                       % judgement parameter
% Picmean = mean2 ( pic ) ;                           % mean of original pic
Residual = zeros ( Times ) ;  Residual ( 1 ) = sum ( abs ( R - SysMatrix * Display ) ) ;      % used as stop condition
figure  % hold residual graph

while ( IterativeTime <= Times && Residual ( IterativeTime ) > 1e-6 )            % end condition of loop
%              disp ( IterativeTime ) ;
             Err = R - SysMatrix * Display ;
             Display = Display + Lamda * SysMatrix' * ( Err ./ ( sum( SysMatrix , 2 ) + EPS ) ) ./ ( sum ( SysMatrix ) + EPS )' ; 
             Display ( Display < 0 ) = 0 ;       % non-negation constraint
%     Dismean = mean ( Display ) ;
    ME ( IterativeTime ) = norm ( Display - picvector ) /  ( Norm * height * width ) ;      % compute error
    IterativeTime = IterativeTime + 1 ;
    
    Residual ( IterativeTime ) = norm ( R - SysMatrix * Display ) / ( Norm * height * width ) ;        % used as stop condition
    plot ( 2 : IterativeTime , Residual ( 2  : IterativeTime ) ) ;
    ylim ( [ 0 , ( 10 * Residual ( IterativeTime ) ) ] ) ;
    drawnow ; 
    
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
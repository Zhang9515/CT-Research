% 2018/03/05 by ZXZ
tic
clear all;
% close all;
%%
% parameter define
thetaint = 3 ;                                                                 % theta unit 
thetaRange = thetaint : thetaint : 180 ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

pic = phantom ( 257 ) ;
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.2 ;

t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , Ltheta * Lt , 1 ) ;       %  to be consistent with sysmatix

%% compute system matrix

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;

picvector = reshape ( flipud( pic )' , height * width , 1  ) ;  % original image
R = SysMatrix * picvector ;        % generate projection with system matrix
% load A.mat
% SysMatrix = A ; 
%% iterative section

Display = zeros ( height * width , 1 ) ;
% Picvector = reshape ( Display) ; 

times = 1 ; 
Stop = 200 ;

% compute tk

% x = ones ( Lt * Ltheta , 1 ) ; 
% M = SysMatrix * SysMatrix' ;
% err = 1 ;
% tmax  = 0 ;
% while ( err >= 1e-6 )
%     
%     tmax0 = tmax ; 
%     x1m = M * x ;
%     x1 = x1m / max ( x1m ) ;   
%     x = x1 ; 
%     tmax = 1 / max ( x1m ) ;  
%     disp ( tmax )
%     err = norm ( abs ( tmax - tmax0 ) ) ;
%     
% end

tmax =  0.0010536 ;
tk = 0.95 * tmax ; 

while ( times <= Stop )
    
    disp ( times ) ;
    % steepest gradient algorithm
    
    BackPro = SysMatrix * Display ; 
    Residual = BackPro - R ; 
    Gradient = 2 * SysMatrix' * Residual ; 
    Display = Display - tk * Gradient ;  
    
    % steepest gradient algorithm ( reduce Matrix-Vector computation )
    
%     if ( mod ( times , 10 ) == 1 )
%         Residual = 2 * SysMatrix' * ( R - SysMatrix * Display ) ;
%     else 
%         Residual = Residual - alpha * SysMatrix' * AG * 2 ;      
%     end
%     AG = SysMatrix * Residual ; 
%     alpha = Residual' * Residual / ( AG' * AG * 2 ) ;
%     Display = Display + alpha * Residual ; 
    
    ME ( times ) = sum ( abs (  Display - picvector ) ) /  ( height * width ) ;    % compute error
    
    times = times + 1 ; 

end
Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;
figure , imshow ( Display , [ 0  1 ] ) ;                     % display results

figure, plot ( 1 : times - 1 , ME( 1  : times - 1 ) ) ;                          % display error graph

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;

toc
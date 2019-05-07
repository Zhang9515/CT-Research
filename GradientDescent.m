% 2018/01/15 by ZXZ
% Gradient descent algorithm
tic
clear all;
% close all;
%%
% parameter define
thetaint = 1 ;                                                                 % theta unit 
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

r0 = SysMatrix' * ( R - SysMatrix * Display ) ;   % initial residual for conjugated gradient algorithm ( CG )
d0 = r0 ;              % initial search direction for conjugated gradient algorithm ( CG )

while ( times <= Stop )
    
    disp ( times ) ;
    % steepest gradient algorithm
    
%     BackPro = SysMatrix * Display ; 
%     Residual = BackPro - R ; 
%     Gradient = 2 * SysMatrix' * Residual ; 
%     AG = SysMatrix * Gradient ;
%     Lamda = ( Gradient' * Gradient ) / ( AG' * AG * 2 ) ; 
%     Display = Display - Lamda * Gradient ;  
    
    % steepest gradient algorithm ( reduce Matrix-Vector computation )
    
%     if ( mod ( times , 10 ) == 1 )
%         Residual = 2 * SysMatrix' * ( R - SysMatrix * Display ) ;
%     else 
%         Residual = Residual - alpha * SysMatrix' * AG * 2 ;      
%     end
%     AG = SysMatrix * Residual ; 
%     alpha = Residual' * Residual / ( AG' * AG * 2 ) ;
%     Display = Display + alpha * Residual ; 
 
    % conjugated gradient algorithm ( CG )
    
    if  ( times == 1 )
            d = d0 ;  
            r_previous = r0 ; 
    else
            r_previous = r_next ; 
    end
    AG = SysMatrix * d ; 
    alpha = ( d' * r_previous ) / ( AG' * AG ) ;
    Display = Display + alpha * d ;
    if ( mod ( times , 10 ) == 1 )
        r_next = SysMatrix' * ( R - SysMatrix * Display ) ;
    else 
        r_next = r_previous - alpha * SysMatrix' * AG ;      
    end 
    beta = ( r_next' * r_next ) / ( r_previous' * r_previous ) ;  
    d = r_next + beta * d ; 
   
    MSE ( times ) = sum ( sum ( ( Display - picvector ).^2 ) ) /  ( height * width ) ;    % compute error
    
    times = times + 1 ; 
    
end
Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;
figure , imshow ( Display , [ 0  1 ] ) ;                     % display results

figure, plot ( 1 : times - 1 , MSE( 1  : times - 1 ) ) ;                          % display error graph

figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
title ( ' grey distrubition ' ) ;
axis ( [ 0 256 0 1 ] ) ;

toc
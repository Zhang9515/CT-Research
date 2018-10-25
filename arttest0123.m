clear all;
% close all;
%%
% parameter define
thetaint = 1 ;                                                                 % theta unit 
thetaRange = thetaint : thetaint : 360 ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

pic = phantom ( 257 ) ;
Size = [ 257 , 257 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 


t_int = 1 ;
t_range =  -184 : t_int : 184 ;
Lt = length ( t_range ) ;

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
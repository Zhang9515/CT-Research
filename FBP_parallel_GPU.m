% modified on 2017/03/16 ZXZ
% modified on 2017/04/7 
tic
%% clear all;
% close all;
%%
% parameter define

% pic = single(StandardPhantom ( "modified_shepp_logan" , 257 ) );
load 'E:\ZXZ\Data\trial2D'
pic = single(trial2D) ; 

Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.1 ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

thetaint = deg2rad(0.5) ;     % theta unit 
Maxtheta = deg2rad(360) ;      
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range, radian
Ltheta = length ( thetaRange ) ; 

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection

%% GPU accelerated Siddon projection

picvector = single( reshape (pic, height * width , 1) ) ;
% t_range = single(t_range') ; thetaRange = single( thetaRange' ) ;
R = ProjectionParallel_2D( picvector , height , width , Size , single(thetaRange') , single(t_range') ) ;     % store parallel beam projection
% R = ProjectionParallel_2D( picvector , height , width , Size , thetaRange , t_range ) ;
% R = reshape( single(R) , Lt , Ltheta ) ;
% figure,imshow(R,[])

Display = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
% Display = reshape( Display , Lt , Ltheta ) ;
Display = reshape( Display , height , width ) ;
Display = Display' ;
figure,imshow(Display , [0 0.5])
          
%%   display

figure,plot( 1 : height , Display ( 257 , : ) , 1 : height , pic ( 257 , : ) ) ;
title ( ' grey distrubition ' ) ;
axis ( [ 0 512 0 2 ] ) ;

% aver=sum(sum(p))/(size(p,1)*size(p,2));
% pd=double(p);
% d=(sum(sum((pd-f).^2))/sum(sum((pd-aver).^2)))^0.5;
toc ;
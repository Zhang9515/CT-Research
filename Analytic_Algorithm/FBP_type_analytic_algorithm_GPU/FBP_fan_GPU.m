%2019/09/14 ZXZ
% equal spaced GPU based
% reference : Algorithms for Reconstruction with Nondiffracting Sources
%              |                     coordinate system in matlab different from the one in the theory
%              |
%              |
% -----------|------------> j
%              |
%              |
%             \/  i
%  **** RScan between source and center point play a great role in suppressing unintersted area *****
%  
%  % we define the actual length of the pic as ( x = 60 mm y = 60 mm )
%   note that resolution is different from length
tic 
clear all;
Size = [ 60 , 60 ] ;                                  % actual range
% pic = single(phantom( 512 )) ;           % original picture : number means the number of pixel in the range 
load('E:\ZXZ\Data\trial2D.mat')
pic = trial2D ;
Displaywindow = [0 0.5] ;

BetaScanInt = deg2rad(5) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project

MaxP = RPic * ( 1 + 0.1 )  ;                                           
PInt = Resolution ;                      %   interval of S ( interval on the detect plain ), empircally pixel-detector ratio is related to size of image
Pdomain = - MaxP : PInt : MaxP ;                          % detective range
LP = length ( Pdomain ) ;

Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  

Ratio = 4 ;                                                           % should be smaller than 8
RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

% R = zeros ( LBeta ,  LP ) ;   % create space to store fan projection
%%  GPU accelerated projection
picvector = Img2vec_Mat2Cpp2D( pic ) ;
R = ProjectionFan_2D ( picvector, height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% R = reshape( R , LP , LBeta )' ;
% figure,imshow(R, [])

Display = FBPfan( single(R) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;
% Display = BackprojectionFan2D( single(R) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;

%% Display 
Display = Vec2img_Cpp2Mat2D( Display , height , width ) ;
figure,imshow( Display , Displaywindow)
% ProjectionDisplay = flipud ( Projection' ) ;
% figure , imshow( ProjectionDisplay , [ 0 , 1 ] ) ; 
% title ( ' Reconstructed image ' ) ;
% pic = phantom( height , width ) ;
% mid_index = round(height/2) ;
% figure , plot ( 1 : size ( pic , 1 ) , ProjectionDisplay ( mid_index , : ) , 1 : height , pic ( mid_index , : ) ) ;
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 height 0 1 ] ) ;
%%   frequency display
% Projectionfft = fft2 ( Projection ) ;
% Projectionffts = fftshift ( Projectionfft ) ;
% Fm = abs ( Projectionffts ) ;  
% figure, imshow ( log ( 1  + Fm ) , [ ] )

%% evaluation function
% aver=sum ( sum ( pic ) ) / ( size ( pic , 1 ) * size ( pic , 2 ) ) ;
% pd=double ( pic ) ;
% d= ( sum ( sum ( ( pd - Projection ).^2 ) ) / sum ( sum ( ( pd - aver ).^2 ) ) ) ^0.5 ;
toc
%2019/03/21
% implement primary dictionary learning with OMP & K_SVD
clear ;
%% image input
load 'E:\ZXZ\Data\trial2D'
pic = single(trial2D) ; 
% figure,imshow(pic,[0 0.5])
% pic = single(StandardPhantom( 256 )) ;
%% parameter define

imgsize = size( pic ) ; 
[ height , width ] = size ( pic ) ;              % store the size of picture

Size = [ 60 , 60 ] ;                                  % actual range

Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

thetaint = deg2rad(1) ;                                                                 % theta unit 
Maxtheta = deg2rad(360) ;
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range
thetaRange = single( thetaRange' ) ; 
Ltheta = length ( thetaRange ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.5 ;
t_range = -tmax : t_int : tmax ;
t_range = single( t_range' ) ; 
Lt = length ( t_range ) ;

% R = zeros ( Lt ,  Ltheta ) ;   
picvector = reshape (pic, height * width , 1) ;
R = ProjectionParallel_2D( picvector , height , width , Size , thetaRange , t_range ) ;     % store parallel beam projection
R = reshape( R , Lt , Ltheta ) ;
figure,imshow(R,[])
%% patch operation

imgsize = size( pic ) ;
patchsize = [5 , 5] ; 
slidestep = [3 , 3] ;
patchset = ExtractPatch2D ( pic , patchsize , slidestep) ; 
Image2D = PatchSynthesis ( patchset , patchsize , slidestep, imgsize ) ;
figure,imshow(Image2D,[0 0.5])
% Alpha = omp( Dictionary , Xintm , Dictionary' * Dictionary , sparsity ) ;



 





















% modified on 2017/03/16 ZXZ
% modified on 2017/04/7 

clear ;
tic
% close all;
%%
% parameter define

% pic = single(StandardPhantom ( 'modified shepp-logan' , 512 ) );
% server2 path
load 'E:\ZXZ\Data\trial2D'
% lab computer path
% load 'G:\CTcode\Data\trial2D'
pic = single(trial2D) ; 
pic = single(phantom(512)) ;

Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = Resolution * 1.1 ;
t_range =  -tmax : t_int : tmax ;
% tmax = 364 ; 
% t_int = 1 ;
% t_range =  (-tmax : t_int : tmax) * Resolution ;
Lt = length ( t_range ) ;

thetaint = deg2rad(1) ;     % theta unit 
Maxtheta = deg2rad(360) ;      
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range, radian
Ltheta = length ( thetaRange ) ; 

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection

%% GPU accelerated Siddon projection

picvector = Img2vec_Mat2Cpp2D ( pic ) ;
% t_range = single(t_range') ; thetaRange = single( thetaRange' ) ;
R = ProjectionParallel_2D( picvector , height , width , Size , thetaRange' , t_range' ) ;     % store parallel beam projection

% for hamtest = 1 : 20
%% hamming
% R = reshape( R , Lt , Ltheta ) ;
% Hamming = zeros ( ( Lt * 2 - 1 )  , 1 ) ; 
% HamRadius = 6 ;
% Hamming ( Lt - HamRadius : Lt + HamRadius ) = hamming ( 2 * HamRadius +1 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% Rcov = zeros( Lt , Ltheta ) ;
% for i = 1 : Ltheta
%         cov = conv ( R ( : , i ) , Hamming ) ;                              % convolution with filter
%         Rcov ( : , i ) = cov ( Lt : 2 * Lt - 1 ) / Hammingsum ;
% end
% Rcov = reshape( Rcov , Lt * Ltheta , 1 ) ;
% figure,imshow(R,[])
    % R1 = radon(trial2D,rad2deg(thetaRange)) * Resolution ;
    % R1 = reshape( R1 ,  Lt * Ltheta ,1 ) ; 
%%
Display = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
% disp(HamRadius) ;
disp( psnr(Display,double(picvector),1) );
% end% hamtest

Display = Vec2img_Cpp2Mat2D( Display , height , width ) ;
% figure,imshow( Display , [0 0.5])
figure,imshow(abs(Display - double(pic)),[])
% %%   display
% 
% figure,plot( 1 : height , Display ( 257 , : ) , 1 : height , pic ( 257 , : ) ) ;
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 512 0 1.2 ] ) ;
% 
% % aver=sum(sum(p))/(size(p,1)*size(p,2));
% % pd=double(p);
% % d=(sum(sum((pd-f).^2))/sum(sum((pd-aver).^2)))^0.5;
toc ;
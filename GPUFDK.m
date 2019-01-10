% 2019/01/10 by ZXZ
% GPUFDK, all codes are extracted from threeDprojection.m, just to make the
% file clear and simple
% ( t , s , z )

%                 t /\     /|  Z
%                    |   /
%   S               | /
% <--------------|---------------
%                  / |
%                /   |
%              /     |
%
%              |                     coordinate system in matlab different from the one in the theory
%              |
%              |
% -----------|------------> j
%              |
%              |
%             \/  i

tic 
clear;
Size = [ 512 * 0.7480 ; 512 * 0.7480 ; 211 * 1 ] ;     % actual range 60

load('E:\ZXZ\Data\ThoraxHD.mat')
pic = single( ThoraxHD ) ;

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture

Resolution = max ( Size ) / t_length ;                           

Distance = 700 ;              % distance between source and center point

PInt = 1 ;                                    %interval of P ( 0.1 exact )

MaxP = ( t_length / 2 * Distance ) / ( Distance - t_length / 2 ) ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

XigamaInt = 1 ;                                      % interval of Xigama ( 0.1 exact )

MaxXigama = ( z_length / 2 * Distance ) / ( Distance - t_length * sqrt(2) / 2 ) ;       % computed by rule of similar triangle
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = single(Xigamadomain');
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

% FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = deg2rad(1) ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = deg2rad(360) ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

%% GPU accelerate projection

picvector = reshape (pic, t_length * s_length * z_length, 1);
clear pic;
% size of R is ( LP * LXigama * LBeta, 1 )
R = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance);

%% GPU accelerate the preweight & filteration & backprojection

z_length = 211;
R= single(R) ;
Display = FDK ( R , Xigamadomain , Pdomain , BetaScanRange , Distance, Size, t_length, s_length, z_length) ;
Display = reshape ( Display , t_length , s_length , z_length ) ;
figure,imshow3Dfull(Display, [0 0.5])

toc


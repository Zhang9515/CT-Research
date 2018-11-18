% 2018/11/14
clear all
tic
Size = [ 64 , 64 , 64 ] ;     % actual range 60

pic = phantom3d ( 'Modified Shepp-Logan' , 64 ) ;     % original picture  
% load ( ' D:\TestCpp\CT\Data\FDK\SFBPimage3.mat ' ) ;       % head ct data
% pic = double ( SFBPimage ) ;
% load ( ' D:\TestCpp\CT\Data\FDK\Chest256_160.mat ' ) ;       % chest ct data
% pic = double ( Chest ) ;

pic = single(pic);

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
%t_length = 64 ; s_length = 64 ; z_length = 64 ;              % store the size of picture
PicSize = [ t_length , s_length , z_length ] ;

Resolution = Size ./ PicSize ;                           
Rplane = max ( Size ) / 2 ;                    % radius of project in the plane

% PInt = 0.4 ;                                    %interval of P ( 0.1 exact )
PInt = 1 ;
MaxP = Rplane * 1.1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

% XigamaInt = 0.4 ;                                      % interval of Xigama ( 0.1 exact )
XigamaInt = 1 ; 
MaxXigama = Rplane * 1.1 ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = single(Xigamadomain');
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

DisRatio = 4 ; 
Distance = 103.9; % Rpic * DisRatio ;                     % distance between source and center point

FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = 1.5 ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = 360 ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

% CovarianceMatrix = eye( LBeta * LXigama * LP ) ;

%% 

picvector = reshape (pic, t_length * s_length * z_length, 1) ;
R = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;

% R = reshape ( R , LP , LXigama, LBeta ) ;
% figure,imshow3Dfull(R,[])



toc
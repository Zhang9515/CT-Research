%% SIRT3D 18/03/15 ZXZ

tic 
clear all ;
Size = [ 64 , 64 , 64 ] ;     % actual range 60
pic = phantom3d ( 'Modified Shepp-Logan' , 64 ) ;     % original picture  
% load ( ' D:\TestCpp\CT\Data\FDK\SFBPimage3.mat ' ) ;       % head ct data
% pic = double ( SFBPimage ) ;
% load ( ' D:\TestCpp\CT\Data\FDK\Chest256_160.mat ' ) ;       % chest ct data
% pic = double ( Chest ) ;

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
%t_length = 64 ; s_length = 64 ; z_length = 64 ;              % store the size of picture
PicSize = [ t_length , s_length , z_length ] ;

Resolution = Size ./ PicSize ;                           
Rplane = max ( Size ) / 2 ;                    % radius of project in the plane

% PInt = 0.4 ;                                    %interval of P ( 0.1 exact )
PInt = 1 ;
MaxP = Rplane * 1.1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
LP = length ( Pdomain ) ;

% XigamaInt = 0.4 ;                                      % interval of Xigama ( 0.1 exact )
XigamaInt = 1 ; 
MaxXigama = Rplane * 1.1 ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

DisRatio = 4 ; 
Distance = 103.9; % Rpic * DisRatio ;                     % distance between source and center point

FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = 1.5 ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = 360 ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix

SysMatrix = GenSysMatCircular3D ( PicSize , Center_t , Center_s , Center_z , Pdomain  , Xigamadomain , BetaScanRange , Distance , Resolution ) ;
for num = 1 : prod ( PicSize )
            zindex = floor ( ( num - 1 ) / ( PicSize(1) * PicSize(2) ) ) + 1 ; 
            sindex = floor ( ( num - ( zindex - 1 ) * PicSize(1) * PicSize(2) - 1 ) / PicSize (1) ) + 1 ; 
            tindex = mod ( num -1 , PicSize (1) ) + 1 ; 
            picvector ( num ) = pic ( tindex , sindex , zindex ) ;
end
picvector = picvector' ; % column vector
R = SysMatrix * picvector;        % generate projection with system matrix
Norm = sum ( R ) ;
disp ( ' System Matrix complete ! ' )
% load A.mat
% SysMatrix = A ; 
%R_superviser = reshape(R, LP, LXigama, LBeta);
%% iterative
% S.Kaczmarz Method

Times = 50 ;
IterativeTime = 1  ;      % times to iterative
Lamda = 1 ;                       % relaxtion factor  SIRT
Display = zeros ( prod ( PicSize ) , 1 ) ;          % store the reconstruction
% Display = ones ( 1  , height * width) ;          % store the reconstruction for MART

ME = zeros( 1 ,  Times ) ;                                       % judgement parameter
% Picmean = mean2 ( pic ) ;                           % mean of original pic
Residual = zeros ( Times ) ;  Residual ( 1 ) = sum ( abs ( R - SysMatrix * Display ) ) ;      % used as stop condition
figure  % hold residual graph

while ( IterativeTime <= Times )            % end condition of loop
%              disp ( IterativeTime ) ;
             Err = R - SysMatrix * Display ;
             Display = Display + Lamda * SysMatrix' * ( Err ./ ( sum( SysMatrix , 2 ) + 1e-9 ) ) ./ ( sum ( SysMatrix ) + 1e-9 )' ; 
             Display ( Display < 0 ) = 0 ;       % non-negation constraint
%     Dismean = mean ( Display ) ;
%     ME ( IterativeTime ) = sum ( abs ( Display - picvector ) ) ./  ( sum * prod ( PicSize ) ) ;      % compute error
    IterativeTime = IterativeTime + 1 ;
    
    Residual ( IterativeTime ) = sum ( abs ( R - SysMatrix * Display ) ) / ( Norm * prod ( PicSize ) ) ;        % used as stop condition
    plot ( 2 : IterativeTime , Residual ( 2  : IterativeTime ) ) ;
    ylim ( [ 0 , ( 10 * Residual ( IterativeTime ) ) ] ) ;
    drawnow ; 
    
end
Display = reshape ( Display , PicSize(1) , PicSize(2) , PicSize(3) ) ;

figure , imshow3Dfull ( Display , [ 0 1 ] ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;

 toc
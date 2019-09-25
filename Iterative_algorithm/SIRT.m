%% SIRT 17/12/16 ZXZ
% update at 19/05/04
% actually what i do is SIRT , which is different from SART in the term of projection line. 
tic
clear all;
% close all;
%%
% lab computer
% load 'G:\CTcode\Data\trial2D'
% server2 path
% parameter define
load ..\..\Data\Adaptive_patchsize_selection\trial2D

% %% fan beam model
% displaywindow = [0 0.5] ;
% 
% BetaScanInt = deg2rad(5) ;             % scanning internal              
% MaxBeta = deg2rad(360) ; 
% BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
% LBeta = length ( BetaScanRange ) ; 
% 
% pic = trial2D ; 
% % clear trial2D
% % pic = phantom(512) ;
% Size = [ 60 , 60 ] ;                                  % actual range
% 
% [ height , width ] = size ( pic ) ;              % store the size of picture
% Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
% RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project
% 
% Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  
% 
% MaxP = RPic * ( 1 + 0.1 )  ;                                           
% PInt = Resolution ;                      %   interval of S ( interval on the detect plain ), empircally pixel-detector ratio is related to size of image
% Pdomain = - MaxP : PInt : MaxP ;                          % detective range
% LP = length ( Pdomain ) ;
% 
% Ratio = 4 ;                                                           % should be smaller than 8
% RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 
% 
% R = zeros ( LP ,  LBeta ) ;   % create space to store fan projection
%% parallel beam model
pic = single(trial2D) ; 

Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

% tmax = round ( Rpic * 1.1 ) ;
% t_int = 0.05 ;
% t_range =  -tmax : t_int : tmax ;
tmax = RPic * ( 1 + 0.1 ) ; 
t_int = Resolution ;
t_range =  (-tmax : t_int : tmax) ;
Lt = length ( t_range ) ;

thetaint = deg2rad(5) ;     % theta unit 
Maxtheta = deg2rad(360) ;      
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range, radian
Ltheta = length ( thetaRange ) ; 

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection

%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
% load ..\..\Data\SysMatrix_parallel_512
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
% tic
% R = SysMatrix * double(picvector) ;        % generate projection with system matrix
% toc
% R = reshape( R , Lt , Ltheta ) ;
% figure,imshow(R,[])
% Norm_pic = norm ( picvector ) ;
% % load A.mat
% % SysMatrix = A ; 
%% GPU-based projection 
% picvector = Img2vec_Mat2Cpp2D( pic ) ;
% R = ProjectionFan_2D ( picvector, height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% tic
R = ProjectionParallel_2D( picvector , height , width , Size , thetaRange' , t_range' ) ; 
% toc
% R2 = reshape( R2 , Lt , Ltheta ) ;
% figure,imshow(R2,[])

%% iterative
% S.Kaczmarz Method
Threshold = 1e-6 ;
Times = 100 ;
IterativeTime = 1  ;      % times to iterative
Lamda = 1 ;                       % relaxtion factor  SIRT
Display = zeros ( height * width , 1 ) ;          % store the reconstruction
MinLim = 0 ; MaxLim = 1;     % range constraint of solution
% Display = ones ( 1  , height * width) ;          % store the reconstruction for MART
LDisplay = length ( Display ) ;
LR = length(R) ;
rmse = zeros( 1 ,  Times ) ;                                       % judgement parameter
% Display_previous = FBPfan( single(R) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;
% Err = R - ProjectionFan_2D ( single(Display_previous), height, width, Size, BetaScanRange', Pdomain', RScan ) ;
Display_previous = FBPparallel( single(R) , single(thetaRange') , single(t_range') , Size , height ,width ) ;
Err = R - ProjectionParallel_2D( single(Display_previous) , height , width , Size , thetaRange' , t_range' ) ;

local_e = LocalError( Display , Display_previous ) ;
figure  % hold residual graph

% compute the row/column sum of system matrix
Rowsum = sum( SysMatrix , 2 ) ;
Colsum = sum ( SysMatrix ) ;
% Rowsum = ProjectionFan_2D ( single(ones(LDisplay,1)), height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% Colsum = BackprojectionFan2D( single(ones(LR,1)) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ;


while ( IterativeTime <= Times && local_e > Threshold )            % end condition of loop
%              
           
             Display = Display_previous + Lamda * (SysMatrix' * ( Err ./ ( Rowsum + eps ) )) ./ ( Colsum + eps )' ; 
%              tem = Err ./ ( Rowsum + eps ) ;
%              Display = Display_previous + Lamda * BackprojectionFan2D( single(tem) , single(BetaScanRange') , single(Pdomain') , Size , height, width, RScan ) ./ ( Colsum + eps ) ; 
             Display ( Display < MinLim ) = MinLim ;       Display ( Display > MaxLim ) = MaxLim ;   % non-negation constraint
%     Dismean = mean ( Display ) ;
            rrmse ( IterativeTime ) = rRMSE( Display , picvector ) ;      % compute error
            local_e = LocalError( Display , Display_previous ) ;
            disp ( ['IterativeTime: ', num2str(IterativeTime), ';   |    rRMSE: ', num2str( rrmse(IterativeTime )),';   |    local_e: ', num2str(local_e) ]) ;
            IterativeTime = IterativeTime + 1 ;
             
%             Display_med = Vec2img_Cpp2Mat2D( Display , height , width ) ;
%             imshow ( Display_med , [ 0  0.5 ] ) ;                     % display results
%             drawnow;
            Display_previous = Display ;
%             Err = R - SysMatrix * Display ;

%     Rmed = ProjectionParallel_2D( single(Display) , height , width , Size , deg2rad(thetaRange') , t_range' ) ;
%     Residual ( IterativeTime ) = norm ( R - Rmed ) / ( Norm * height * width ) ;  
    
%     plot ( 2 : IterativeTime , Residual ( 2  : IterativeTime ) ) ;
%     ylim ( [ 0 , ( 10 * Residual ( IterativeTime ) ) ] ) ;
%     drawnow ; 
    
end

Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;

figure , imshow ( Display , [ 0  0.5 ] ) ;                     % display results

% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;


 toc
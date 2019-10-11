% 2018/10/19 
% a new fast splitting algorithm to solve the Weighted Split Bregman minimization problem 
% in the backward step of an accelerated Forward-Backward algorithm
% AFB: accelerated Forward-Backward algorithm
  
tic
clear all;
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D
displaywindow = [0 0.5] ;
% close all;
%% parallel beam
% parameter define
% thetaint = 1 ;                                                                 % theta unit 
% thetaRange = thetaint : thetaint : 180 ;                                % radon scanning range
% Ltheta = length ( thetaRange ) ; 
% 
% pic = phantom ( 128 ) ;
% % pic = dicomread ( 'ActualCTImage.dcm');
% % winL = 0 ;    winH = 4095 ;           %  set window width
% % pic = winL + double ( pic ) / ( winH - winL ) ;        
% 
% Size = [ 60 , 60 ] ;                                  % actual range
% 
% [ height , width ] = size ( pic ) ;              % store the size of picture
% Resolution = max ( Size ) / height ;   % define the resolution of pic
% Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
% Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 
% 
% tmax = round ( Rpic * 1.1 ) ;
% t_int = 0.2 ;
% t_range =  -tmax : t_int : tmax ;
% Lt = length ( t_range ) ;
% 
% R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection
%% fan beam
% parameter define
BetaScanInt = deg2rad(5) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

pic = trial2D ; 
% clear trial2D
% pic = phantom(512) ;
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project

Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  

MaxP = RPic * ( 1 + 0.1 )  ;                                           
PInt = Resolution ;                      %   interval of S ( interval on the detect plain ), empircally pixel-detector ratio is related to size of image
Pdomain = - MaxP : PInt : MaxP ;                          % detective range
LP = length ( Pdomain ) ;

Ratio = 4 ;                                                           % should be smaller than 8
RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

R = zeros ( LP ,  LBeta ) ;   % create space to store fan projection
%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix
load ..\..\..\Data\SysMatrix_fan_512
% SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
picvector = Img2vec_Mat2Cpp2D( pic ) ;  % original image
R = SysMatrix * double(picvector) ;        % generate projection with system matrix
Display = SysMatrix' * R;  % initiate the iterate output
miu = 0.1 ;  % weight of the data fidelity term
alpha = 0.1 ;      % this parameter is the balance weight between TV and prior image
% Matrix = @(v) miu * SysMatrix' * (SysMatrix * v) + v * ( 1 - alpha );
% [matrixEigen , iterate] = powermethod ( Matrix , randn( height * width , 1 ) ) ; 

%% accelerated forward-backward iterate
matrixEigen = 494.1557 * miu +2*(1-alpha) ;    % largest eigenvalue of a specific system matrix
step = 0.9 / matrixEigen ;
a = 2; % used as parameter for modified FISTA
t_previous = ( 0 + a + 1) / a;
MinLim = 0 ; MaxLim = 1 ; 
% N_outloop = 20 ; 
% i_outloop = 1;
Display_previous = Display ; 
Display_prior = double(picvector) ;
gradientMatrix_x = gradient2Dmatrix_x(height,width);
gradientMatrix_y = gradient2Dmatrix_y(height,width);

divergence_matrix = divergenceMatrix2D(height,width);
delta = norm ( divergence_matrix , inf ) ;
lamda = 0.1 / ( step * delta ) ;        % ALM paramter, influence the convergence of the inner loop in the split bregman

for i_outloop  = 1 : 10
    grad_G = miu * SysMatrix' * ( R - SysMatrix * Display_previous ) + 2 * (1 - alpha) * ( Display_previous - Display_prior )  ;
    V = Display_previous + step * grad_G ;
%     [ delta_x , delta_y , delta_wx , delta_wy , wtv ] = WTV ( Display , height , width ) ;

    U = FWSB ( gradientMatrix_x , gradientMatrix_y , V , alpha , lamda , step, MinLim , MaxLim ) ;
    
    t = ( i_outloop + a +1) / a;
    if i_outloop ==1
        Display = U ;
    else
        tao = ( t_previous - 1) / t; 
        Display = U_previous + tao * ( U - U_previous ) ;
    end
    local_error = LocalError ( Display , Display_previous ) ;
    loss = alpha * (norm(gradientMatrix_x * Display,1) + norm(gradientMatrix_y * Display,1)) + ( 1 - alpha ) * (norm(Display-Display_prior))...
                 + 0.5 * miu * norm( SysMatrix * Display - R  ) ;   
    disp( ['i_outloop: ' , num2str(i_outloop) , ' local_error: ', num2str(local_error),' global_error: ', num2str(LocalError ( Display , picvector )),...
        ' loss: ' , num2str(loss)] ) ;
      
    U_previous = U;
    t_previous = t;
    Display_previous = Display ;
    test = Vec2img_Cpp2Mat2D( Display , height , width ) ;
    figure,imshow(test,displaywindow ) ;
    drawnow ; 
end

Display = Vec2img_Cpp2Mat2D( Display , height , width ) ;
figure,imshow(Display,displaywindow)

toc
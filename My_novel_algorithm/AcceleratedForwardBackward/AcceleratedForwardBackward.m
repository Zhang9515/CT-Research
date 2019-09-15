% 2018/10/19 
% a new fast splitting algorithm to solve the Weighted Split Bregman minimization problem 
% in the backward step of an accelerated Forward-Backward algorithm
% AFB: accelerated Forward-Backward algorithm
  
tic
clear all;
% close all;
%%
% parameter define
thetaint = 1 ;                                                                 % theta unit 
thetaRange = thetaint : thetaint : 180 ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

pic = phantom ( 128 ) ;
% pic = dicomread ( 'ActualCTImage.dcm');
% winL = 0 ;    winH = 4095 ;           %  set window width
% pic = winL + double ( pic ) / ( winH - winL ) ;        

Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.2 ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection
%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
picvector = reshape ( flipud( pic )' , height * width , 1  ) ;  % original image
R = SysMatrix * picvector ;        % generate projection with system matrix
Display = SysMatrix' * R;  % initiate the iterate output
a = SysMatrix' * SysMatrix ; 

% [lamdaEigen , iterate] = powermethod ( SysMatrix' * SysMatrix , randn( height * width , 1 ) ) ; 

lamdaEigen = 11469.3363107980 ;    % largest eigenvalue of a specific system matrix
beta = 0.1 / lamdaEigen ;
a = 2; % used to parameter for modified FISTA
t_previous = ( 0 + a + 1) / a;
r = 10e-4;

% N_outloop = 20 ; 
i_outloop = 1;
Display_previous = Display ; 
while (1)
    disp( [ 'i_outloop:' , num2str(i_outloop) ] ) ;
    lamda = r * norm( Display ,1 ) ;
    V = Display + beta * SysMatrix' * ( R - SysMatrix * Display );
    [ delta_x , delta_y , delta_wx , delta_wy , wtv ] = WTV ( Display , height , width ) ;
    delta = norm ( delta_wx' * delta_wx + delta_wy' * delta_wy , inf ) ;
    theta = 0.5 / ( beta * delta ) ;
    U = FWSB ( delta_wx , delta_wy , V , lamda , theta , beta ) ;
    
    t = ( i_outloop + a +1) / a;
    if i_outloop ==1
        Display = U ;
    else
        alpha = ( t_previous - 1) / t; 
        Display = U_previous + alpha * ( U - U_previous ) ;
    end
    
    if ( norm ( Display_previous - Display , 2) <= 10e-6 * norm( Display_previous , 2 ) )
          break;
    end
    disp( norm ( Display_previous - Display , 2) / norm( Display_previous , 2 ) )
    disp( norm ( picvector - Display , 2) / norm( picvector , 2 ) )
    U_previous = U;
    t_previous = t;
    Display_previous = Display ;
    i_outloop = i_outloop + 1 ; 
end

Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;
figure,imshow(Display,[])

toc
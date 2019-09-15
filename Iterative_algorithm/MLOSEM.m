% MLEM & OSEM 2018/01/24 by ZXZ
tic
clear all;
% close all;
%%
% parameter define
thetaint = 0.5 ;                                                                 % theta unit 
thetaRange = thetaint : thetaint : 180 ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

% pic = phantom ( 257 ) ;
pico = dicomread ( 'ActualCTImage.dcm');
winL = 0 ;    winH = 4095 ;           %  set window width
pico = winL + double ( pico ) / ( winH - winL ) ;        

pic ( 2 : 513 , 2 : 513 ) = pico ; clear pico 
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = Resolution ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;
%% compute system matrix

SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;

picvector = reshape ( flipud( pic )' , height * width , 1  ) ;  % original image
R = SysMatrix * picvector ;        % generate projection with system matrix
Norm = sum ( R ) ;
%% iterative section  MLEM

% Display = ones ( height * width , 1 ) ;     % For EM never all is zero on the denominator 
% % Picvector = reshape ( Display) ; 
% 
% times = 1 ; 
% Stop = 400 ;
% 
% while ( times <= Stop )
%     
%     disp ( times ) ;
% 
%      AG = SysMatrix * Display ; 
%      Scurrent = spdiags ( Display ./ sum ( SysMatrix )' , 0 , height * width , height * width ) ;
%      Lamda = spdiags ( 1 ./ ( AG + 1e-10 ) , 0 , Ltheta * Lt , Ltheta * Lt ) ;
%      Display = Display + Scurrent * ( SysMatrix' * ( Lamda * ( R - AG ) ) ) ;
%      
%     MSE ( times ) = sum ( sum ( ( Display - picvector ).^2 ) ) /  ( height * width ) ;    % compute error
%     
%     times = times + 1 ; 
%     
% end
%% iterative section  OSEM

Display = ones ( height * width , 1 ) ;     % For EM never all is zero on the denominator 
Subset = 10 ;                        % when Subset =1 , it is MLEM
LSub = Ltheta * Lt / Subset ;

times = 1 ; 
Stop = 100 ;
% Residual = zeros ( Stop ) ;  Residual ( 1 ) = sum ( abs ( R - SysMatrix * Display ) ) / ( height * width ) ) ;      % used as stop condition
Condition = zeros ( Stop ) ; Condition ( 1 ) = 1 ; 
ME = zeros ( Stop ) ;  
figure     % hold residual graph

while ( times <= Stop && Condition ( times ) > 1e-6 )
    
%     disp ( times ) ;
    PreDisplay = Display ; 
    for num = 1 : Subset     
        
        SubSysMatrix = SysMatrix ( ( num - 1 ) * LSub + 1 : num * LSub , : ) ;     % split system matrix into smaller subsets
        SubR = R ( ( num - 1 ) * LSub + 1 : num * LSub ) ; 
        AG = SubSysMatrix * Display ; 
        Scurrent = spdiags ( Display ./ ( sum ( SubSysMatrix )' + 1e-10 ) , 0 , height * width , height * width ) ;      % sparse diagnal matrix
        Lamda  = spdiags ( 1 ./ ( AG + 1e-10 ) , 0 , LSub , LSub ) ;                                      % sparse diagnal matrix
        Display = Display + Scurrent * ( SubSysMatrix' * ( Lamda * ( SubR - AG ) ) ) ;
        
    end
     
    ME ( times ) = sum ( abs ( Display - picvector ) ) ./  ( Norm * ( height * width ) ) ;    % compute normalized error
    
    times = times + 1 ; 
%     Residual ( times ) = sum ( abs ( R - SysMatrix * Display ) / ( height * width ) ) ;        % used as stop condition
    Condition ( times ) = sum ( abs ( Display - PreDisplay ) / ( height * width ) ) ;      % Covergence criterion
    
    plot ( 2 : times , Condition ( 2  : times ) ) ;
    ylim ( [ 0 , ( 10 * Condition ( times ) ) ] ) ;
    drawnow ; 
    
end
close

%% Display
Display = reshape ( Display , width , height ) ;
Display =  flipud ( Display' ) ;
figure , imshow ( Display , [ 0  0.5 ] ) ;                     % display results

figure, plot ( 1 : times - 1 , ME( 1  : times - 1 ) ) ;                          % display error graph
title ( ' error graph ' ) ;

% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 256 ) , 1 : size ( pic , 1 ) , pic ( : , 256 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 ( height - 1 ) 0 1 ] ) ;

toc






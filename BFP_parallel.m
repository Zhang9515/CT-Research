% 2017/07/18 ZXZ
% first : BFP in parallel projection
tic
clear all;
% close all;
%%
% parameter define
thetaint = 0.5 ;                                                                 % theta unit 
thetaRange = 90 + thetaint : thetaint : 270 ;                                % radon scanning range
Ltheta = length ( thetaRange ) ; 

pic = StandardPhantom ( 513 ) ;
% pic = phantom ( 513 ) ;
% pic = flipud ( pic ) ;
Size = [ 60 , 60 ] ;                                  % actual range
X1 = -Size ( 1 ) / 2 ; X2 = Size ( 1 ) / 2 ;    Y1 = -Size ( 1 ) / 2 ; Y2 = Size ( 1 ) / 2 ;  % periphral points

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = Size ( 1 ) / height ;   % define the resolution of pic
Center_x = Size ( 1 ) / 2 ; Center_y = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;              % 
t_int = 0.117 ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

R = zeros ( Ltheta , Lt ) ;   % create space to store fan projection
%% formualr
%         A    a     b    x0    y0    phi
%        ---------------------------------
shep = [  2   .69   .92    0     0     0                                                      % in ground coordinate 
        -.98 .6624 .8740   0  -.0184   0
        -.02 .3100 .1100   .22    0    72
        -.02 .4100 .1600  -.22    0     108
         .01 .2100 .2500   0    .35    0
         .01 .0460 .0460   0    .1     0
         .01 .0460 .0460   0   -.1     0
         .01 .0460 .0230 -.08  -.605   0 
         .01 .0230 .0230   0   -.606   0
         .01 .0460 .0230   .06  -.605   90   ] ;

% shep = [  1   .69   .92    0     0     0                                                      % in ground coordinate
%         -.8 .6624 .8740   0  -.0184   0
%         -.2 .3100 .1100   .22    0    72
%         -.2 .4100 .1600  -.22    0     108
%          .1 .2100 .2500   0    .35    0
%          .1 .0460 .0460   0    .1     0
%          .1 .0460 .0460   0   -.1     0
%          .1 .0460 .0230 -.08  -.605   0 
%          .1 .0230 .0230   0   -.606   0
%          .1 .0460 .0230   .06  -.605   90   ] ;
     
proportion = max ( Size ) / 2 ;                            % actual scale
R = reshape ( R , 1 , [] ) ; 

parfor num = 1 : Ltheta * Lt
    
        tindex = floor ( ( num - 1 ) / Ltheta ) + 1 ;
        thetaindex = mod ( num - 1 ,  Ltheta ) + 1 ;
        thetaRadian =  thetaRange ( thetaindex ) * pi / 180 ;       % angle to radian            
        t = t_range ( tindex ) ;      
                    for n = 1 : 10
                        x0 = shep ( n , 4 ) * proportion ; y0 = shep ( n , 5 ) * proportion ;              % oval information in ground coordinate
                        s = sqrt ( x0^2 + y0^2 ) ;
                        ggama = atan2 ( y0 , x0 ) ;
                        alpaha = shep ( n , 6 ) * pi / 180 ;
                        A = shep ( n , 2 ) * proportion ; B = shep ( n , 3 ) * proportion  ;    
                        Rou = shep ( n , 1 ) ;                   
                        ts = t - s * cos ( ggama - thetaRadian ) ;
                        thetas = thetaRadian - alpaha ; 
                        a2 = ( A * cos ( thetas ) )^2 + ( B * sin ( thetas ) )^2 ;
                        if ( abs ( ts ) <= sqrt ( a2 ) )                                             % decide whether cross the oval
                                R ( num ) = R ( num ) + 2 * Rou * A * B * sqrt ( a2 - ts^2 ) / a2 ;
                        end                             
                    end   
end
R = reshape ( R , Ltheta , Lt ) ;

%% differentiate 

RD =  zeros ( Ltheta ,  Lt ) ;

for i  = 1 : Ltheta
        for j  = 1 : Lt
                if ( j ~= Lt )
                    RD ( i , j ) = ( R ( i ,  j ) - R ( i , j + 1 ) ) / t_int ; 
%                     RD ( i , j ) = ( R ( i ,  j + 1) - R ( i , j ) ) / t_int ; 
                else
                    RD ( i , j ) = ( R ( i ,  Lt - 1 ) - R ( i , Lt ) ) / t_int ; 
%                     RD ( i , j ) = ( R ( i ,  Lt ) - R ( i , Lt  - 1 ) ) / t_int ; 
                end                    
        end  
end   

%% backprojection on PI-line

Rate = 0.4 ; 
Lpad = round ( Rate * ( width - 1 )  / 2 ) ;                                % length of pad
PHeight = height ; PWidth = width + 2 * Lpad ;
ProjectionShift = zeros ( Ltheta , PHeight , PWidth ) ;
ProjectionShift = reshape ( ProjectionShift , 1 , [] ) ;
Resolution2 = max ( Size ) / width ;
ResolutionPI = Resolution2 ;           % resolution of PI line
LPI = Size ( 1 ) + 2 * Lpad * Resolution2 ;                                % length of PI line , a little bit larger than ROI
XPI1 = -LPI / 2 ;    XPI2 = LPI / 2 ;
PCenter_x = LPI / 2 ;
G = zeros ( PHeight , PWidth ) ;                        % store back-projection

parfor num = 1 : Ltheta * PHeight * PWidth                       
            
            y = floor ( ( num - 1 ) / ( PWidth * Ltheta ) ) + 1 ;                % piexl number
            x = floor ( ( num - ( y - 1 ) * PWidth * Ltheta - 1 )  / Ltheta ) + 1 ; 
            k = mod ( num - 1 , Ltheta ) + 1 ;   
            
            thetaRadian = thetaRange ( k ) * pi / 180 ;
            Point_y = ( y - 0.5 ) * Resolution2 ; Point_x = ( x - 0.5 ) * Resolution2 ;            % center point of the pixel
            
            t =  ( Point_x - PCenter_x ) * cos ( thetaRadian ) + ( Point_y - Center_y ) * sin ( thetaRadian ) ;
            t1index = floor ( ( t + tmax ) / t_int ) + 1 ;
            t2index = t1index + 1 ;
            if ( t >= min ( t_range ) && t < max ( t_range ) ) 
                ProjectionShift ( num ) = ( ( t_range ( t2index ) - t ) * RD ( k , t1index ) + ( t - t_range ( t1index ) ) * RD ( k , t2index ) ) / t_int * thetaint * pi / 180 ;  % interpolation
            end
end

for i = 1 : PHeight
    for j = 1 : PWidth
            G ( i , j ) = sum ( ProjectionShift ( 1 + ( i - 1 ) * PWidth * Ltheta + ( j - 1 ) * Ltheta : Ltheta + ( i - 1 ) * PWidth * Ltheta + ( j - 1 ) * Ltheta ) ) ;
    end    
end       

%% Hilbert filtration
% PI-line is horizontal projection line, Hilbert function is less accurate in the peripheral area. 
% Solution 1: increase the sampling rate ; Solution 2: add zeros pad at two ends ( when LPI is about 1.5 width , R reach minimum) 

Dpi =  R ( end , : ) ;            Dpi = fliplr ( Dpi ) ;     %  projection along PI line ( theta= 270 , direction of t and i is the same )
% Dpi =  sum ( pic' ) * Resolution ;

Display = zeros ( height , width ) ;   % store reconstruction result
Display = reshape ( Display , 1 , [] ) ;

% PIwidth = round ( Size ( 1 ) / ResolutionPI ) ; 
% Hilb = zeros ( height , width ) ;                          % store Hibert filtration    
XPIRange =  XPI1 + ( linspace( 1 , PWidth , PWidth ) - 0.5 ) * ResolutionPI ; 

parfor num = 1 : height * width   

            Windex = floor ( ( num - 1 ) / height ) + 1 ;                          % pixel number
            Hindex = mod ( num - 1  , height ) + 1 ; 
            Hilb = [] ;
            XPI = X1 + ( Windex - 0.5 ) * Resolution ; 
%             for PIindex = 1 : PIwidth 
%                     xpi = X1 + ( PIindex - 0.5 ) * ResolutionPI ; 
%                     if ( xpi ~= XPI )
%                         Hilb ( PIindex ) = sqrt ( ( XPI2 - xpi ) * ( xpi - XPI1 ) ) / ( XPI - xpi ) * G ( Hindex , Windex ) * ResolutionPI ;    
%                     end
%             end
            HindexG = round ( Hindex * Resolution / ResolutionPI ) ;
            WindexG = round ( Windex * Resolution / ResolutionPI ) ;
            Hilb = imag( hilbert ( sqrt ( ( XPI2 - XPIRange ) .* ( XPIRange - XPI1 ) ) .* G ( HindexG , : ) ) ) ;         % hilbert tranform 
%            Hilb = imag( hilbert ( G ( HindexG , : ) ) ) ; 
           
            Y = Y1 + ( Hindex  - 0.5 ) * Resolution ; 
            tindex = round ( ( Y - t_range ( 1 ) ) / t_int ) + 1 ;
            tindex1 = floor ( ( Y - t_range ( 1 ) ) / t_int ) + 1 ; tindex2 = tindex1 + 1 ;
            t1 =  t_range ( tindex1 ) ; t2 =  t_range ( tindex2 ) ; 
            DPI = ( ( t2 - Y ) * Dpi ( tindex1 ) + ( Y - t1 ) * Dpi ( tindex2 ) ) / ( t_int ) ; 
%             DPI = Dpi ( Hindex ) ;
%             Display ( num ) = (  sum ( Hilb )  +   2 * pi * DPI ) / ( 2 * pi ^ 2 * sqrt ( ( XPI2 - XPI ) * ( XPI - XPI1 ) ) ) ; 
%             Display ( num ) = (  Hilb  ( WindexG + Lpad ) * pi  +   2 * pi * DPI ) / ( 2 * pi ^ 2 * sqrt ( ( XPI2 - XPI ) * ( XPI - XPI1 ) ) ) ; 
            Display ( num ) =  (  Hilb  ( WindexG + Lpad ) ) / ( 2 * pi ) ; 

end
clear G ;
Display = ( reshape ( Display , height , width ) ) ;

% figure , imshow ( flipud ( Display ) , [ 0 , 1 ] ) ; 
figure , imshow ( flipud ( Display ) , [ 1,1.05 ] ) ; 
Display =  flipud ( Display ) ;
figure,plot( 1 : size ( pic , 1 ) , Display ( : , 257 ) , 1 : size ( pic , 1 ) , pic ( : , 257 ) ) ;
title ( ' grey distrubition ' ) ;
axis ( [ 0 512 1 1.05 ] ) ;
% R computation

% Rindex = sum ( sum ( abs ( Display - pic ) ) ) / sum ( sum ( abs ( pic ) ) ) ;

toc ;

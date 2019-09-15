% 2017/08/04 ZXZ
% second : BFP in fanbeam projection
tic
% clear all;
% close all;
%%
% parameter define
Betaint = 1 ;                                                                 % theta unit 
BetaScanRange = 0 : Betaint : 360 - Betaint ;                                % radon scanning range
LBeta = length ( BetaScanRange ) ; 

% pic = StandardPhantom ( 513 ) ;
pic = phantom ( 513 ) ;
pic = flipud ( pic ) ;
Size = [ 60 , 60 ] ;                                  % actual range
X1 = -Size ( 1 ) / 2 ; X2 = Size ( 1 ) / 2 ;    Y1 = -Size ( 1 ) / 2 ; Y2 = Size ( 1 ) / 2 ;  % periphral points

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = Size ( 1 ) / height ;   % define the resolution of pic
Center_x = Size ( 1 ) / 2 ; Center_y = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

MaxS = Rpic * ( 1 + 0.1 )  ;                                           
SInt = 0.1 ;                      %   interval of S ( interval on the detect plain )
Sdomain = - MaxS : SInt : MaxS ;                          % detective range
LS = length ( Sdomain ) ;

Ratio = 4 ;                                                           % should be smaller than 8
RScan = Rpic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

R = zeros ( LBeta , LS ) ;   % create space to store fan projection
%% formualr
%         A    a     b    x0    y0    phi
%        ---------------------------------
% shep = [  2   .69   .92    0     0     0                                                      % in ground coordinate
%         -.98 .6624 .8740   0  -.0184   0
%         -.02 .3100 .1100   .22    0    72
%         -.02 .4100 .1600  -.22    0     108
%          .01 .2100 .2500   0    .35    0
%          .01 .0460 .0460   0    .1     0
%          .01 .0460 .0460   0   -.1     0
%          .01 .0460 .0230 -.08  -.605   0 
%          .01 .0230 .0230   0   -.606   0
%          .01 .0460 .0230   .06  -.605   90   ] ;

shep = [  1   .69   .92    0     0     0                                                      % in ground coordinate
        -.8 .6624 .8740   0  -.0184   0
        -.2 .3100 .1100   .22    0    72
        -.2 .4100 .1600  -.22    0     108
         .1 .2100 .2500   0    .35    0
         .1 .0460 .0460   0    .1     0
         .1 .0460 .0460   0   -.1     0
         .1 .0460 .0230 -.08  -.605   0 
         .1 .0230 .0230   0   -.606   0
         .1 .0460 .0230   .06  -.605   90   ] ;
     
proportion = max ( Size ) / 2 ;                            % actual scale
R = reshape ( R , 1 , [] ) ; 

parfor num = 1 : LBeta * LS 
    
        Betaindex = mod ( num - 1 ,  LBeta ) + 1 ;
        Sindex = floor ( ( num - 1 ) / LBeta ) + 1 ;
        
        beta =  BetaScanRange ( Betaindex ) - 90 ;                                              % in my view beta is the angle between projection lay and x-axis
        betaRadian = beta * pi / 180 ;                    % angle to radian
        %define the source in matlab coordinate
        source_x = Center_x + RScan * cos ( betaRadian ) ;  source_y = Center_y + RScan * sin ( betaRadian ) ;         % in matlab coordinate
                   
                    gamaRadian = atan ( Sdomain ( Sindex ) / RScan ) ;
                    t = RScan * sin ( gamaRadian ) ;
                    thetaRadian = gamaRadian + betaRadian ;
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
R = reshape ( R , LBeta , LS ) ;

%% differentiate
RD =  zeros ( LBeta ,  LS ) ;

for i  = 1 : LBeta
        if ( i > 181 )
            signal = 1 ; 
        else 
            signal = -1 ; 
        end
        for j  = 1 : LS
                if ( j ~= LS )
                    RD ( i , j ) = signal * ( R ( i ,  j ) - R ( i , j + 1 ) ) / SInt ; 
                else
                    RD ( i , j ) = signal * ( R ( i ,  LS - 1 ) - R ( i , LS ) ) / SInt ; 
                end                    
        end  
end   

%% backprojection on PI-line
G = zeros ( size ( pic ) ) ;                        % store back-projection

G = reshape ( G , 1 , [] ) ;
Resolution2 = max ( Size ) / height ;

parfor num = 1 : height * width                       
            
            Windex = floor ( ( num - 1 ) / height ) + 1 ;                          % pixel number
            Hindex = mod ( num - 1  , height ) + 1 ;             
            Point_x = X1 + ( Windex - 0.5 ) * Resolution2 ;  Point_y = Y1 + ( Hindex - 0.5 ) * Resolution2 ;            % center point of the pixel
            r = [ Point_x , Point_y ] ;
            Lamda = asin ( Point_y / RScan ) * 180 / pi ; 
            if ( Lamda >= 0 )
                    lamdaRange = floor ( Lamda ) : ceil ( 180 - Lamda ) ;                                          % PI-line scanning angular range, angular interval = 1
            else
                    lamdaRange = floor ( 180 - Lamda  ) : ceil ( 360 + Lamda ) ;
            end
            Llamda = length ( lamdaRange ) ;
            g = zeros ( Llamda , 1 ) ; 
            for lamdaindex = 1 : Llamda
                lamdaRadian = lamdaRange ( lamdaindex ) * pi / 180 ;
                u = - Point_x * sin ( lamdaRadian ) + Point_y * cos ( lamdaRadian ) ;                      % rotation coordinate 
                w = Point_x * cos ( lamdaRadian ) + Point_y * sin ( lamdaRadian ) - RScan ;           
                ud = RScan * u / w ;                                                                                     % projection coordinate , ud-xis on the detector-line is converted to u
                r0 = RScan .* [ ( cos ( lamdaRadian ) ) , ( sin ( lamdaRadian ) ) ] ;      % source in fixed coordinate
                Dr0 = RScan .* [ ( -sin ( lamdaRadian ) ) , ( cos ( lamdaRadian ) ) ] ;   % derivation of r0
                Distance = norm ( r - r0 , 2 ) ;                                                                      % distance between point and source 
                eU = [ ( -sin ( lamdaRadian ) ) , ( cos ( lamdaRadian ) ) ] ;                                  % unit vector of axis-U
                Orientation = ( r - r0 ) ./ Distance ;
                A = norm ( [ ud , RScan ] , 2 ) ;
                udIndex1 = floor ( ( ud - Sdomain ( 1 ) ) / SInt ) + 1 ; udIndex2 = udIndex1 + 1 ;
                if ( ud >= min ( Sdomain ) && ud < max ( Sdomain ) )
                        ud1 = Sdomain ( udIndex1 ) ;  ud2 = Sdomain ( udIndex2 ) ;
                        if ( lamdaRange ( lamdaindex ) >= 360 )
                            Lamdaindex = int16 ( lamdaRange ( lamdaindex ) - 360 ) + 1 ;
                        else 
                            Lamdaindex = int16 ( lamdaRange ( lamdaindex ) ) + 1 ; 
                        end
                        P = ( ( ud2 - ud ) * R ( Lamdaindex , udIndex1 ) + ( ud - ud1 ) * R ( Lamdaindex , udIndex2 ) ) / SInt ; 
                        DP = ( ( ud2 - ud ) * RD ( Lamdaindex , udIndex1 ) + ( ud - ud1 ) * RD ( Lamdaindex , udIndex2 ) ) / SInt ; 
                end               
                g ( lamdaindex ) = ( - ( Dr0 * Orientation' ) * P + ( Dr0 * eU' ) * A * DP ) / Distance ^ 2 * Betaint * pi / 180 ;
                if ( lamdaindex == 1 )
                        g1 = P / Distance ; 
                elseif ( lamdaindex == Llamda )
                        g2 = P / Distance ; 
                end
            end
            G ( num ) = sum ( g ) + g2 - g1 ;
end
G = reshape ( G , height , width ) ;

%% Hilbert filtration
% PI-line is horizontal projection line
% load G ; load R ; 
Dpi = zeros ( height , 1 ) ;                                          %  projection along PI line ( theta= 90 , direction of t and i is the same )
parfor num = 1 : height
        Point_y = Y1 + ( num - 0.5 ) * Resolution2 ; 
        Lamda = asin ( Point_y / RScan ) ; 
        s = - RScan * tan ( Lamda ) ;
        sIndex = floor ( ( s - Sdomain ( 1 ) ) / SInt ) + 1 ;
        LamdaAng = Lamda * 180 / pi ;
        if ( LamdaAng >= -0.5 )
            LamdaIndex  = round ( LamdaAng ) + 1 ;
        else 
            LamdaIndex = round ( 360 + LamdaAng ) + 1 ;
        end
        Dpi ( num ) = R ( LamdaIndex , sIndex ) ;
end

LPI = Size ( 1 ) ;                                % length of PI line , a little bit larger than ROI
XPI1 = -LPI / 2 ;    XPI2 = LPI / 2 ;
Display = zeros ( height , width ) ;   % store reconstruction result
Display = reshape ( Display , 1 , [] ) ;
ResolutionPI = Resolution ; 
PIwidth = round ( Size ( 1 ) / ResolutionPI ) ; 
% Hilb = zeros ( height , width ) ;                          % store Hibert filtration    
XPIRange =  X1 + ( linspace( 1 , width , width ) - 0.5 ) * Resolution ; 

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
            Hilb = imag( hilbert ( sqrt ( ( XPI2 - XPIRange ) .* ( XPIRange - XPI1 ) ) .* G ( Hindex , : ) ) ) ;         % hilbert tranform 
           
%             Display ( num ) = (  sum ( Hilb )  +   2 * pi * DPI ) / ( 2 * pi ^ 2 * sqrt ( ( XPI2 - XPI ) * ( XPI - XPI1 ) ) ) ; 
                Display ( num ) = (  Hilb  ( Windex ) * pi  +   2 * pi * Dpi ( Hindex ) ) / ( 2 * pi ^ 2 * sqrt ( ( XPI2 - XPI ) * ( XPI - XPI1 ) ) ) ; 

end
Display = ( reshape ( Display , height , width ) ) ;

figure , imshow ( flipud ( Display ) , [ 0 , 1 ] ) ; 
figure,plot( 1 : size ( pic , 1 ) , Display ( : , 257 ) , 1 : size ( pic , 1 ) , pic ( : , 257 ) ) ;
title ( ' grey distrubition ' ) ;
axis ( [ 0 512 0 1 ] ) ;
toc ;

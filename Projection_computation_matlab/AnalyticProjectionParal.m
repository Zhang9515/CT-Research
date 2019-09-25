% 2018/01/16 by ZXZ
% formular parallel projection , mode 1 :  Low Contrast ; mode 2 :  High Contrast
% output : row index is theta and column index is t 

function [ R ] = AnalyticProjectionParal ( Size , thetaRange , t_range , mode )

switch mode 
    
    case 1    %   This is the default head phantom, taken from AK Jain, 439.  ( Low Contrast )
%                      A    a     b    x0    y0    phi
%                ---------------------------------
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
        
    case 2   %   This is the default head phantom, taken from AK Jain, 439.  ( High Contrast )
%                      A    a     b    x0    y0    phi
%                ---------------------------------
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
end
     
proportion = max ( Size ) / 2 ;                            % actual scale
Ltheta = length ( thetaRange ) ;
Lt = length ( t_range ) ;
R = zeros ( Lt ,  Ltheta ) ;   % create space to store parallel projection
R = reshape ( R , 1 , [] ) ; 

parfor num = 1 : Ltheta * Lt

        thetaindex = mod ( num - 1 ,  Ltheta ) + 1 ;
        tindex = floor ( ( num - 1 ) / Ltheta ) + 1 ;
        
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



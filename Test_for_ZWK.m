% test for ZWK 2017/03/16 
% ZXZ
%3D projection Feldkamp
% ( t , s , z )
%    phantom(64) , maxP= 0.75 , pint = 1 , maxxigama = 0.75 , xigamaint = 1 , times = 616.527 s
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
clear all;
BetaScanInt = 1 ;             % scanning internal              
MaxBeta = 360 ; 
pic = phantom3d ( 64 ) ;     % original picture    
BetaScanRange = 0 : BetaScanInt : MaxBeta ;     % scanning range , angle between SO and aixs S
[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
PInt = 1 ;                                    %interval of S
MaxP = max ( size ( pic ) ) * 1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
XigamaInt = 1 ;                                      %interval of S
MaxXigama = max ( size ( pic ) ) * 1 ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Center_t = t_length / 2 ;  Center_s = s_length / 2 ;   Center_z =  z_length / 2 ;          % define the center 
DisRatio = 1 ; 
Distance = max ( size ( pic ) ) * DisRatio ;                     % distance between source and center point

%% direct fan : equal spaced
% 257.438 s  phantom(64) , maxP= 0.75 , pint = 1 
R = zeros ( t_length ,  s_length , z_length ) ;   % create space to store fan projection
distanceInt = 1 ;      %define the interval on the rayline
for i = 1 :  length ( BetaScanRange ) 
    beta =  BetaScanRange ( i ) ;
    betaRadian = beta * pi / 180 ;                    % angle to radian
    source_t = Center_t - Distance * sin ( betaRadian ) ;  source_s = Center_s + Distance * cos ( betaRadian ) ;  source_z  = Center_z ;     %define the source in matlab coordinate
    for j = 1 : length ( Xigamadomain ) 
        gama = atan ( Xigamadomain ( j ) / Distance ) ;          % radian angle in s-z coordinate plane
        Distance_shift = Distance / cos ( gama ) ;                    % length of DO'
        for k = 1 : length ( Pdomain )
              theta = atan ( Pdomain ( k ) / Distance_shift ) ;                                    % radian angle in s'-t coordinate plane 

                %   Siddon algrithim
                Smax = 2 * Distance ; 
                DetectPoint_tend = Center_t + Smax * sin ( theta ) * cos ( betaRadian ) - ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * sin ( betaRadian ) ; 
                DetectPoint_send = Center_s + Smax * sin ( theta ) * sin ( betaRadian ) + ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * cos ( betaRadian ) ;
                DetectPoint_zend = Center_z + Smax * cos ( theta ) * sin ( gama ) ;
                if ( DetectPoint_tend > source_t )                                                         % define projection line direction and range of i and j
                     t_signal = 1 ;   
                     t_range = ceil ( source_t ) : floor ( DetectPoint_tend ) ;
                else 
                    t_signal = -1 ; 
                    t_range = ceil ( DetectPoint_tend ) : floor ( source_t ) ;
                end
                if ( DetectPoint_send > source_s )
                    s_signal = 1 ; 
                    s_range = ceil ( source_s ) : floor ( DetectPoint_send ) ;
                else
                    s_signal = -1 ; 
                    s_range = ceil ( DetectPoint_send ) : floor ( source_s ) ;
                end     
                if ( DetectPoint_zend > source_z )
                    z_signal = 1 ; 
                    z_range = ceil ( source_z ) : floor ( DetectPoint_zend ) ;
                else
                    z_signal = -1 ; 
                    z_range = ceil ( DetectPoint_zend ) : floor ( source_z ) ;
                end
                ProjectionLine = zeros ( 3 , length ( t_range ) + length ( s_range ) + length ( z_range ) ) ;
                ProjectionLine ( 1 , 1 : length ( t_range ) ) = abs ( ( t_range - source_t ) / ( DetectPoint_tend - source_t + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , 1 : length ( t_range ) ) = t_range ;
                ProjectionLine ( 3 , 1 : length ( t_range ) ) = 1 ;                              % 1 represents t
                ProjectionLine ( 1 , ( length ( t_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) ) ) = abs ( ( s_range - source_s ) / ( DetectPoint_send - source_s + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( length ( t_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) ) ) = s_range ;
                ProjectionLine ( 3 , ( length ( t_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) ) ) = 2 ;                              % 2 represents s
                ProjectionLine ( 1 , ( length ( t_range ) + length ( s_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) + length ( z_range ) ) ) = abs ( ( z_range - source_z ) / ( DetectPoint_zend - source_z + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( length ( t_range ) + length ( s_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) + length ( z_range ) ) ) = z_range ;
                ProjectionLine ( 3 , ( length ( t_range ) + length ( s_range ) + 1 ) : ( length ( s_range ) + length ( t_range ) + length ( z_range ) ) ) = 3 ;                              % 3 represents z
                ProjectionLine = sortrows ( ProjectionLine' )' ;
                DetectPoint_t = round ( source_t ) ; DetectPoint_s = round ( source_s ) ; DetectPoint_z = round ( source_z ) ;                   % start point of the line        
                for n = 1 : length ( ProjectionLine )
                    if ( ProjectionLine ( 3 , n ) == 1 )
                        if  ( t_signal == 1 )
                            DetectPoint_t = ProjectionLine ( 2 , n ) + 1 ; 
                        elseif ( t_signal == -1 )
                            DetectPoint_t = ProjectionLine ( 2 , n ) ; 
                        end
                    elseif ( ProjectionLine ( 3 , n ) == 2 )
                        if  ( s_signal == 1 )
                            DetectPoint_s = ProjectionLine ( 2 , n ) + 1 ; 
                        elseif ( s_signal == -1 )
                            DetectPoint_s = ProjectionLine ( 2 , n ) ; 
                        end
                    elseif ( ProjectionLine ( 3 , n ) == 3 )
                        if  ( z_signal == 1 )
                            DetectPoint_z = ProjectionLine ( 2 , n ) + 1 ; 
                        elseif ( z_signal == -1 )
                            DetectPoint_z = ProjectionLine ( 2 , n ) ; 
                        end
                    end
                    if ( DetectPoint_s > 0 && DetectPoint_s <= s_length && DetectPoint_t > 0 && DetectPoint_t <= t_length && DetectPoint_z > 0 && DetectPoint_z <= z_length ) 
                        LengthUnit =  ( ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ) ;  
                        R ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) = R ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + LengthUnit ;    
                    end
                end
        end
    end
end

% figure , imagesc( squeeze ( R ( : , ( length ( Xigamadomain ) + 1 ) / 2  , : ) ) ) ; 
% save D:\TestCpp\R.mat R ;
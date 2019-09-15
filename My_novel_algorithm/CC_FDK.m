%2017/3/13 ZXZ
% CC-FDK
% ( t , s , z )
%    phantom(64) , maxP= 0.75 , pint = 1 , maxxigama = 0.75 , xigamaint = 1 , times = 616.527 s
%                 s /\     /|  Z
%                    |   /
%   t               | /
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
BetaScanInt = deg2rad(1) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
pic = phantom3d ( 128 ) ;     % original picture    
BetaScanRange = 1 : BetaScanInt : MaxBeta ;     % scanning range , angle between SO and aixs S

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture

PInt = 0.5 ;                                    %interval of P
MaxP = max ( size ( pic ) ) * 1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
XigamaInt = 1 ;                                      % interval of Xigama
MaxXigama = max ( size ( pic ) ) * 1 ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Center_t = t_length / 2 ;  Center_s = s_length / 2 ;   Center_z =  z_length / 2 ;          % define the center 
DisRatio = 2 ; 
Distance = max ( size ( pic ) ) * DisRatio ;                     % distance between source and center point

%% direct fan : equal spaced
% 257.438 s  phantom(64) , maxP= 0.75 , pint = 1 
R = zeros ( length ( BetaScanRange ) ,  length ( Xigamadomain ) , length ( Pdomain ) ) ;   % create space to store fan projection
distanceInt = 1 ;      %define the interval on the rayline
for i = 1 :  length ( BetaScanRange ) 
    beta =  BetaScanRange ( i ) ;
    betaRadian = beta * pi / 180 ;                                          % angle to radian
    source_t = Center_t - Distance * sin ( betaRadian ) ;      %define the source in matlab coordinate
    source_s = Center_s + Distance * cos ( betaRadian ) ;  
    source_z  = Center_z ;                                                        
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
                Lt = length ( t_range ) ; Ls = length ( s_range ) ; Lz = length ( z_range ) ; 
                ProjectionLine = zeros ( 3 , Lt + Ls + Lz ) ;
                ProjectionLine ( 1 , 1 : Lt ) = abs ( ( t_range - source_t ) / ( DetectPoint_tend - source_t + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , 1 : Lt ) = t_range ;
                ProjectionLine ( 3 , 1 : Lt ) = 1 ;                              % 1 represents t
                ProjectionLine ( 1 , ( Lt + 1 ) : ( Ls + Lt ) ) = abs ( ( s_range - source_s ) / ( DetectPoint_send - source_s + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( Lt + 1 ) : ( Ls + Lt ) ) = s_range ;
                ProjectionLine ( 3 , ( Lt + 1 ) : ( Ls + Lt ) ) = 2 ;                              % 2 represents s
                ProjectionLine ( 1 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = abs ( ( z_range - source_z ) / ( DetectPoint_zend - source_z + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = z_range ;
                ProjectionLine ( 3 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = 3 ;                              % 3 represents z
                ProjectionLine = sortrows ( ProjectionLine' )' ;
                DetectPoint_t = ceil ( source_t ) ; DetectPoint_s = ceil ( source_s ) ; DetectPoint_z = ceil ( source_z ) ;                   % start point of the line     
                
                for n = 1 : length ( ProjectionLine )
                    
                    if ( ProjectionLine ( 3 , n ) == 1 && Lt ~= 0 )                                                               % define the pixel 
                            DetectPoint_t = ProjectionLine ( 2 , n ) + 0.5 + t_signal * 0.5 ; 
                    elseif ( ProjectionLine ( 3 , n ) == 2 && Ls ~= 0 )
                            DetectPoint_s = ProjectionLine ( 2 , n ) + 0.5 + s_signal * 0.5 ; 
                    elseif ( ProjectionLine ( 3 , n ) == 3 && Lz ~= 0 )
                            DetectPoint_z = ProjectionLine ( 2 , n ) + 0.5 + z_signal * 0.5 ;       
                    end
                   
                    if ( Lt == 0 && source_t == 0 && DetectPoint_s > 0 && DetectPoint_s <= s_length ...                          % boundary condition
                            && DetectPoint_z > 0 && DetectPoint_z <= z_length )   
                        
                        f = pic ( 1 , DetectPoint_s , DetectPoint_z ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , j , k ) = R ( i , j , k ) + f * LengthUnit ;    
                    end   
                    if ( Ls == 0 && source_s == 0 && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
                            && DetectPoint_z > 0 && DetectPoint_z <= z_length )   
                        
                        f = pic ( DetectPoint_t , 1 , DetectPoint_z ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , j , k ) = R ( i , j , k ) + f * LengthUnit ;    
                    end      
                    if ( Lz == 0 && source_z == 0 && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
                            && DetectPoint_s > 0 && DetectPoint_s <= s_length )    
                        
                        f = pic ( DetectPoint_t , DetectPoint_s , 1 ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , j , k ) = R ( i , j , k ) + f * LengthUnit ;    
                    end   
                    
                    if ( DetectPoint_s > 0 && DetectPoint_s <= s_length && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
                            && DetectPoint_z > 0 && DetectPoint_z <= z_length ) 

                         if ( Lt == 0 && source_t == DetectPoint_t )
                             if ( source_t == t_length )
                                 f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t + 1 , DetectPoint_s , DetectPoint_z ) ) / 2 ;
                             end
                         elseif ( Ls == 0 && source_s == DetectPoint_s ) 
                             if ( source_s == s_length )
                                 f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t , DetectPoint_s + 1 , DetectPoint_z ) ) / 2 ;
                             end
                         elseif ( Lz == 0 && source_z == DetectPoint_z ) 
                             if ( source_z == z_length )
                                 f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z + 1 ) ) / 2 ;
                             end
                         else 
                                f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) ;
                         end
                         
                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                         R ( i , j , k ) = R ( i , j , k ) + f * LengthUnit ;    
                    end
                end
               R ( i , j , k ) = R ( i , j , k ) * Distance / sqrt ( Distance^2 + Pdomain ( k )^2 + Xigamadomain ( j )^2 ) ;
        end
    end
end
figure , imshow ( squeeze ( R ( : , ( length ( Xigamadomain ) + 1 ) / 2  , : ) ) , [] ) ; 
% save D:\TestCpp\Rfdk.mat R ;
 
 %% filter 
G = zeros ( length ( Pdomain )  , 1 ) ;                                            % create filter
for i = 1 : length ( Pdomain ) 
    N = i - ( length  ( Pdomain ) + 1 ) / 2 ;
        if N == 0
            G ( i ) = 1 / ( 8 * PInt ^2 ) ;    % radian 
        elseif rem ( N , 2 ) == 0          % even
            G ( i ) = 0 ;
        else                                         % odd
           G ( i ) = - 0.5 / ( pi * Pdomain ( i ) )  ^ 2 ;
        end
end
Rcov = zeros ( length ( BetaScanRange ) ,  length ( Xigamadomain ) , length ( Pdomain ) ) ;
for i = 1 :  length ( BetaScanRange )
    for j = 1 : length ( Xigamadomain )
    
%     Rsq = zeros ( 1 , length ( Pdomain ) ) ;
%             for jj = 1 : length ( Pdomain )
%                 Rsq ( jj ) = R ( i , j , jj ) ; 
%             end 
        
        cov = conv ( squeeze ( R ( i , j , : ) ) , G' ) ;                              % convolution with filter
        Rcov ( i , j , : ) = PInt * cov ( round ( length ( Pdomain ) * 0.5 ) : round ( length ( Pdomain ) * 0.5 ) + length ( Pdomain ) - 1 ) ;
       % continuous domain to discret domain
        % Rcov ( i , : ) = cov ( 1 :  length ( GamaDetectRange )  )  ;
        %Rcov ( i , : ) = cov (  length ( GamaDetectRange ) : 2 * length ( GamaDetectRange ) - 1 )  ;
    end
end
Hamming = zeros ( length ( Pdomain )  , 1 ) ; 
Hamming ( length ( Pdomain ) / 2 - 7 : length ( Pdomain ) / 2 + 7 ) = hamming ( 15 ) ;                   % convolve with hamming window
Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
for i = 1 :  length ( BetaScanRange )
        cov = conv ( squeeze ( R ( i , j , : ) ) , Hamming' ) ;                              % convolution with filter
        Rcov ( i , j , : ) = cov ( round ( length ( Pdomain ) * 0.5 ) : round ( length ( Pdomain ) * 0.5 ) + length ( Pdomain ) - 1 ) / Hammingsum ;
end
figure , imshow ( squeeze ( Rcov ( : , ( length ( Xigamadomain ) + 1 ) / 2 , : ) ) , [] ) ; 

%% reconstruct

ProjectionCC = zeros ( t_length , s_length , z_length ) ;      %reconstruct the pic
for t = 1 : t_length
    for s = 1 : s_length
        for z = 1 : z_length 
             for betas = 1 : length ( BetaScanRange )
                    beta =  BetaScanRange ( betas ) ;
                    betaRadian = beta * pi / 180 ;
                    source_t = Center_t + Distance * cos ( betaRadian + pi ) ;                           % define the source
                    source_s = Center_s + Distance * sin ( betaRadian + pi ) ;  
                    source_z = Center_z ; 
                    
                    image_t = t - Center_t - 0.5 ;  image_s = s - Center_s - 0.5 ; image_z = z - Center_z - 0.5 ;           % image pixel in ground coordinate
                    
                    dect_t = image_t * cos ( betaRadian ) + image_s * sin ( betaRadian ) ;          % in rotate coordinate
                    dect_s = - image_t * sin ( betaRadian ) + image_s * cos ( betaRadian ) ; 
                    dect_z = image_z ;     
                    
                    LengthRatio = Distance / ( Distance - dect_s ) ; 
                    Xigama_domain = dect_z * LengthRatio ; 
                    P_domain = dect_t * LengthRatio ; 
                    Xigama_domainN1 =  floor ( ( Xigama_domain + MaxXigama ) / XigamaInt ) + 1 ; 
                    Xigama_domainN2 = Xigama_domainN1 + 1 ;
                    P_domainN1 =  floor ( ( P_domain + MaxP ) / PInt ) + 1 ; 
                    P_domainN2 = P_domainN1 + 1 ;
                    
                    if  ( P_domain >= -MaxP && P_domain <= MaxP && Xigama_domain >= -MaxXigama && Xigama_domain <= MaxXigama )
                            P_domain1 = Pdomain ( P_domainN1 ) ; P_domain2 = Pdomain ( P_domainN2 ) ;
                            Xigama_domain1 = Xigamadomain ( Xigama_domainN1 ) ; Xigama_domain2 = Xigamadomain ( Xigama_domainN2 ) ;  
                            % bilinear interploration
                            Xig1 = Xigama_domain - Xigama_domain1 ; Xig2 = Xigama_domain2 - Xigama_domain ;
                            P1 = P_domain - P_domain1 ; P2 = P_domain2 - P_domain1 ;
                            ProjectionCC ( t , s , z ) = ProjectionCC ( t , s , z ) + ( Xig2 * P2 * Rcov ( betas , Xigama_domainN1 , P_domainN1 ) ...
                                + Xig1 * P2 * Rcov ( betas , Xigama_domainN2 , P_domainN1 ) + Xig2 * P1 * Rcov ( betas , Xigama_domainN1 , P_domainN2 ) ...
                                + Xig1 * P1 * Rcov ( betas , Xigama_domainN2 , P_domainN2 ) ) * LengthRatio^2 / ( PInt * XigamaInt ) ;
                    end
             end
             ProjectionCC ( t , s , z ) = ProjectionCC ( t , s , z  ) * BetaScanInt * pi / 180 ;                  
        end
    end
end

figure , imshow( squeeze ( ProjectionCC ( : ,  : , Center_z ) ) ) ; 

% figure , plot ( 1 : size ( pic , 1 ) , squeeze ( Projection ( size ( pic , 1 ) / 2 , : , 50 ) ) , 1 : size ( pic , 1 ) , squeeze ( pic ( size ( pic , 1 ) / 2 , : , 50 ) )  ) ;
% title ( ' grey distrubition ' ) ;

%% CC-modified

load D:\TestCpp\CT\Data\FDK\FDK64R1 ;                             % introduce inner trajectory result
DisRatio0 = 1 ;
Distance0 = max ( size ( pic ) ) * DisRatio0 ; 
ProjectionCCFDK = ( ProjectionCC - Projection ) * Distance0 ^ 2 / ( Distance ^ 2 - Distance0 ^ 2 ) + ProjectionCC ;
figure , imshow( squeeze ( ProjectionCCFDK  ( : ,  : , Center_z ) ) ) ; 




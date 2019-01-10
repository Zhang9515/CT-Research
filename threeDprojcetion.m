%2016/12/03 ZXZ
%2017/4/8 modified
% 2018/4/14 GPU
%3D projection Feldkamp
% ( t , s , z )

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
tic 
clear;
Size = [ 512 * 0.7480 ; 512 * 0.7480 ; 211 * 1 ] ;     % actual range 60
% pic = phantom3d ( 'Modified Shepp-Logan' , 65 ) ;     % original picture  
% pic = Diskphantom ( 512 ) ;
load('E:\ZXZ\Data\ThoraxHD.mat')
pic = single( ThoraxHD ) ;
% load ( ' D:\TestCpp\CT\Data\FDK\SFBPimage3.mat ' ) ;       % head ct data
% pic = double ( SFBPimage ) ;
% load ( ' D:\TestCpp\CT\Data\FDK\Chest256_160.mat ' ) ;       % chest ct data
% pic = double ( Chest ) ;
% load ('G:\CTcode\Code\ProjectionCone_3D\phantom512.mat')
% pic = single(pic);

% PicSize = 512 ;
[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture

%  t_length = 512 ; s_length = 512 ; z_length = 200 ;
 
% t_length = PicSize ; s_length = PicSize ; z_length = PicSize ;              % store the size of picture
Resolution = max ( Size ) / t_length ;                           
% Rpic = max ( Size ) * sqrt ( 3 ) / 2 ;                                         % radius of project (51.9615 for size 60)
% Rpic = 500 ;

Distance = 500 ;              % distance between source and center point

PInt = 0.7 ;                                    %interval of P ( 0.1 exact )
% PInt = 0.2 ;
% MaxP = Rplane * 1.1 ;
MaxP = ( t_length / 2 * Distance ) / ( Distance - t_length / 2 ) ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

XigamaInt = 1 ;                                      % interval of Xigama ( 0.1 exact )
% XigamaInt = 0.2; 
% MaxXigama = Rplane * 1.1 ;
MaxXigama = ( z_length / 2 * Distance ) / ( Distance - t_length * sqrt(2) / 2 ) ;       % computed by rule of similar triangle
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = single(Xigamadomain');
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

% DisRatio = 4 ; 
% Distance = 207.6; % Rpic * DisRatio ;        Distance = 500 ;              % distance between source and center point
% Distance = 730 ; 


% FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = deg2rad(1) ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = deg2rad(360) ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

R = zeros ( LP * LXigama * LBeta, 1 ) ;   % create space to store fan projection

%% GPU accelerate projection

picvector = reshape (pic, t_length * s_length * z_length, 1);
clear pic;
R = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance);
% R = reshape ( R , LP , LXigama, LBeta ) ;
% figure,imshow3Dfull(R,[])

% load('E:\ZXZ\Data\matlab.mat')
% RR = zeros( 870 , 256 , 360 ) ;
% for i = 1 :360
%     RR( : ,:,i) = Proj_matrix(:,:,i)';
% end
% R = reshape(RR,LP * LXigama * LBeta, 1) ;


% clear picvector;
%% direct fan : equal spaced
% 257.438 s  phantom(64) , maxP= 0.75 , pint = 1 
% 
% R = reshape ( R , 1 , [] ) ;                 % adjust the R format to parallel compute
% distanceInt = 1 ;      %define the interval on the rayline
% parfor num = 1 :  LBeta * LXigama * LP 
%     
%             Pindex = floor ( ( num - 1 ) / ( LBeta * LXigama ) ) + 1 ; 
%             Xigamaindex = floor ( ( num - ( Pindex - 1 ) * LXigama * LBeta - 1 ) / LBeta ) + 1 ; 
%             Betaindex = mod ( num -1 , LBeta ) + 1 ; 
% 
%             beta =  BetaScanRange ( Betaindex ) ;
%             betaRadian = beta * pi / 180 ;                                          % angle to radian
%             source_t = Center_t - Distance * sin ( betaRadian ) ;      %define the source in matlab coordinate
%             source_s = Center_s + Distance * cos ( betaRadian ) ;  
%             source_z  = Center_z ;                                                        
% 
%                 gama = atan ( Xigamadomain ( Xigamaindex ) / Distance ) ;          % radian angle in s-z coordinate plane
%                 Distance_shift = Distance / cos ( gama ) ;                    % length of DO'
% 
%                       theta = atan ( Pdomain ( Pindex ) / Distance_shift ) ;                                    % radian angle in s'-t coordinate plane 
% 
%                         %   Siddon algorithm
%                         Smax = 2 * Distance ; 
%                         DetectPoint_tend = Center_t + Smax * sin ( theta ) * cos ( betaRadian ) - ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * sin ( betaRadian ) ; 
%                         DetectPoint_send = Center_s + Smax * sin ( theta ) * sin ( betaRadian ) + ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * cos ( betaRadian ) ;
%                         DetectPoint_zend = Center_z + Smax * cos ( theta ) * sin ( gama ) ;
%                         if ( DetectPoint_tend > source_t )                                                         % define projection line direction and range of i and j
%                              t_signal = 1 ;   
%                              t_range = ceil ( source_t / Resolution ) : floor ( DetectPoint_tend / Resolution ) ;
%                         else 
%                             t_signal = -1 ; 
%                             t_range = ceil ( DetectPoint_tend / Resolution ) : floor ( source_t / Resolution ) ;
%                         end
%                         if ( DetectPoint_send > source_s )
%                             s_signal = 1 ; 
%                             s_range = ceil ( source_s / Resolution ) : floor ( DetectPoint_send / Resolution ) ;
%                         else
%                             s_signal = -1 ; 
%                             s_range = ceil ( DetectPoint_send / Resolution ) : floor ( source_s / Resolution ) ;
%                         end     
%                         if ( DetectPoint_zend > source_z )
%                             z_signal = 1 ; 
%                             z_range = ceil ( source_z / Resolution ) : floor ( DetectPoint_zend / Resolution ) ;
%                         else
%                             z_signal = -1 ; 
%                             z_range = ceil ( DetectPoint_zend / Resolution ) : floor ( source_z / Resolution ) ;
%                         end
%                         Lt = length ( t_range ) ; Ls = length ( s_range ) ; Lz = length ( z_range ) ; 
%                         ProjectionLine = zeros ( 3 , Lt + Ls + Lz ) ;
%                         ProjectionLine ( 1 , 1 : Lt ) = abs ( ( t_range * Resolution - source_t ) / ( DetectPoint_tend - source_t + 1e-10 ) ) * Smax ;
%                         ProjectionLine ( 2 , 1 : Lt ) = t_range ;
%                         ProjectionLine ( 3 , 1 : Lt ) = 1 ;                              % 1 represents t
%                         ProjectionLine ( 1 , ( Lt + 1 ) : ( Ls + Lt ) ) = abs ( ( s_range * Resolution - source_s ) / ( DetectPoint_send - source_s + 1e-10 ) ) * Smax ;
%                         ProjectionLine ( 2 , ( Lt + 1 ) : ( Ls + Lt ) ) = s_range ;
%                         ProjectionLine ( 3 , ( Lt + 1 ) : ( Ls + Lt ) ) = 2 ;                              % 2 represents s
%                         ProjectionLine ( 1 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = abs ( ( z_range * Resolution - source_z ) / ( DetectPoint_zend - source_z + 1e-10 ) ) * Smax ;
%                         ProjectionLine ( 2 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = z_range ;
%                         ProjectionLine ( 3 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = 3 ;                              % 3 represents z
%                         ProjectionLine = sortrows ( ProjectionLine' )' ;
%                         DetectPoint_t = ceil ( source_t / Resolution ) ; DetectPoint_s = ceil ( source_s / Resolution ) ; DetectPoint_z = ceil ( source_z / Resolution ) ;                   % start point of the line     
% 
%                         for n = 1 : length ( ProjectionLine )
% 
%                             if ( ProjectionLine ( 3 , n ) == 1 && Lt ~= 0 )                                                               % define the pixel 
%                                     DetectPoint_t = ProjectionLine ( 2 , n ) + 0.5 + t_signal * 0.5 ; 
%                             elseif ( ProjectionLine ( 3 , n ) == 2 && Ls ~= 0 )
%                                     DetectPoint_s = ProjectionLine ( 2 , n ) + 0.5 + s_signal * 0.5 ; 
%                             elseif ( ProjectionLine ( 3 , n ) == 3 && Lz ~= 0 )
%                                     DetectPoint_z = ProjectionLine ( 2 , n ) + 0.5 + z_signal * 0.5 ;       
%                             end
% 
%                             if ( Lt == 0 && source_t == 0 && DetectPoint_s > 0 && DetectPoint_s <= s_length ...                          % boundary condition
%                                     && DetectPoint_z > 0 && DetectPoint_z <= z_length )   
% 
%                                 f = pic ( 1 , DetectPoint_s , DetectPoint_z ) / 2 ;
%                                 LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                                 R ( num ) = R ( num ) + f * LengthUnit ;    
%                             end   
%                             if ( Ls == 0 && source_s == 0 && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
%                                     && DetectPoint_z > 0 && DetectPoint_z <= z_length )   
% 
%                                 f = pic ( DetectPoint_t , 1 , DetectPoint_z ) / 2 ;
%                                 LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                                 R ( num ) = R ( num ) + f * LengthUnit ;    
%                             end      
%                             if ( Lz == 0 && source_z == 0 && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
%                                     && DetectPoint_s > 0 && DetectPoint_s <= s_length )    
% 
%                                 f = pic ( DetectPoint_t , DetectPoint_s , 1 ) / 2 ;
%                                 LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                                 R ( num ) = R ( num ) + f * LengthUnit ;    
%                             end   
% 
%                             if ( DetectPoint_s > 0 && DetectPoint_s <= s_length && DetectPoint_t > 0 && DetectPoint_t <= t_length ...
%                                     && DetectPoint_z > 0 && DetectPoint_z <= z_length ) 
% 
%                                  if ( Lt == 0 && source_t / Resolution == DetectPoint_t )
%                                          if ( source_t == t_length )
%                                              f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
%                                          else 
%                                              f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t + 1 , DetectPoint_s , DetectPoint_z ) ) / 2 ;
%                                          end
%                                  elseif ( Ls == 0 && source_s / Resolution == DetectPoint_s ) 
%                                          if ( source_s == s_length )
%                                              f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
%                                          else 
%                                              f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t , DetectPoint_s + 1 , DetectPoint_z ) ) / 2 ;
%                                          end
%                                  elseif ( Lz == 0 && source_z / Resolution == DetectPoint_z ) 
%                                          if ( source_z == z_length )
%                                              f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) / 2 ;
%                                          else 
%                                              f = ( pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) + pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z + 1 ) ) / 2 ;
%                                          end
%                                  else 
%                                         f = pic ( DetectPoint_t , DetectPoint_s , DetectPoint_z ) ;
%                                  end
% 
%                                  LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                                  R ( num ) = R ( num ) + f * LengthUnit ;    
%                             end
%                         end
%                        R ( num ) = R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% end
% R = reshape ( R , LBeta , LXigama , LP ) ;
% figure , imshow ( squeeze ( R ( : , 331  , : ) ) , [] ) ; 
% save D:\TestCpp\Rfdk.mat R ;

%%  formular projection  P.104
%                    A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
% Shepp =    [   2  .6900  .920  .900      0       0       0      0    
%                  -.98  .6624  .874  .880      0       0       0      0    
%                  -.02  .310  .1100  .220    .22      0      -.25     72    
%                 -.02  .410  .1600  .280   -.22       0     -.25    108     
%                  .02  .0460  .046  .046      0      .1     -.25     0     
%                  .02  .2100  .250  .500      0      .35    -.25     0     
%                  .01  .0460  .046  .046      0     -.1     -.25     0      
%                  .01  .0460  .023  .020   -.08    -.605    -.25     90    
%                  .01  .0230  .023  .020      0   -.606     -.25     90   
%                  .01  .0460  .023  .020    .06   -.605     -.25      0  
%                  .02  .0560  .040  .10    .06    -.105     .625      0
%                 -.02  .0560  .056  .10       0     .105    .625      0 ];
% % 
% % Disk phantom
%  %         A      a     b     c     x0      y0      z0    phi 
% %        -----------------------------------------------------
% % Shepp =    [  1    .7  .7   .06      0       0       0        0      
% %                     1    .7  .7   .06      0       0      .24        0
% %                     1    .7  .7   .06      0       0       -.24        0
% %                     1    .7  .7   .06      0       0       .48        0
% %                     1    .7  .7   .06      0       0       -.48        0
% %                     1    .7  .7   .06      0       0       .72        0
% %                     1    .7  .7   .06      0       0       -.72      0 ] ;
% 
% Proportion = max ( Size ) / 2 ;        % actual length
% R = reshape ( R , 1 , [] ) ;                 % adjust the R format to parallel compute
% 
% parfor num = 1 :  LBeta * LXigama * LP 
%             
%             Pindex = floor ( ( num - 1 ) / ( LBeta * LXigama ) ) + 1 ; 
%             Xigamaindex = floor ( ( num - ( Pindex - 1 ) * LXigama * LBeta - 1 ) / LBeta ) + 1 ; 
%             Betaindex = mod ( num -1 , LBeta ) + 1 ; 
%             
%             
%             beta =  BetaScanRange ( Betaindex ) ;
%             betaRadian = beta * pi / 180 ;                                          % angle to radian
%             source_t = Center_t + Distance * cos ( betaRadian + pi / 2 ) ;      % define the source in ground coordinate
%             source_s = Center_s + Distance * sin ( betaRadian + pi / 2 ) ;    % source locates on the S axis
%             source_z  = Center_z ;                                                        
%                     
%                     gama = atan ( Xigamadomain ( Xigamaindex ) / Distance ) ;          % radian angle in s-z coordinate plane ( is converted compared with book )
%                     r = Xigamadomain ( Xigamaindex ) * Distance / sqrt ( Distance^2 + Xigamadomain ( Xigamaindex )^2 ) ;                % parameter in parallel ray
%                     Distance_shift = Distance / cos ( -gama ) ;                    % length of DO'
%                              
%                              FanA = atan ( Pdomain ( Pindex ) / Distance ) ; 
%                              thetaRadian = betaRadian + FanA ;                                    % radian angle in s'-t coordinate plane 
%                              t = Pdomain ( Pindex ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 ) ;                % parameter in parallel ray
%                              for n = 1 : 12
%                                         
%                                         A = Shepp ( n , 2 ) * Proportion ; B = Shepp ( n , 3 ) * Proportion ; C = Shepp ( n , 4 ) * Proportion ;                    % information of oval
%                                         rou = Shepp ( n , 1 ) ;
%                                         X1 = Shepp ( n , 5 ) * Proportion ; Y1 = Shepp ( n , 6 ) * Proportion ; Z1 = Shepp ( n , 7 ) * Proportion ;
%                                         fai = Shepp ( n , 8 ) * pi / 180 ;
%                                         
%                                         thetas = thetaRadian - fai ;
%                                         ts = t - X1 * cos ( thetaRadian ) - Y1 * sin ( thetaRadian ) ;
%                                         rs = r + ( Y1 * cos ( thetaRadian ) - X1 * sin ( thetaRadian ) ) * sin ( -gama ) - Z1 * cos ( -gama ) ;   
%                                         
%                                         a2 = C^2 * ( ( B * sin ( thetas ) )^2 + ( A * cos ( thetas ) )^2 ) * ( cos ( -gama ) )^2 + ( A * B * sin ( -gama ) )^2 ;
%                                         Sq = a2 - ts^2 * ( ( C * cos ( -gama ) )^2 + ( ( B * cos ( thetas ) )^2 + ( A * sin ( thetas ) )^2 ) * ( sin ( -gama ) )^2 ) ...
%                                                 - ( rs )^2 * ( ( B * sin ( thetas ) )^2 + ( A * cos ( thetas ) )^2 ) * ( 7 + cos ( -4 * gama ) ) / 8 ...
%                                                 - 2 * ts * rs * sin ( -gama ) * cos ( thetas ) * sin ( thetas ) * ( B^2 - A^2 ) ;
%                                         if ( Sq >= 0 )
%                                                 R ( num ) = R ( num ) + 2 * rou * A * B * C * sqrt ( Sq ) / a2 ;
%                                         end                            
%                              end
%                              
% % % half-scan-weighted
% %                              if ( betaRadian <= 2 * ( FanAmax - FanA ) )
% %                                     W = sin ( pi * betaRadian / ( FanAmax - FanA ) / 4 )^2 ;
% %                              elseif ( betaRadian > 2 * ( FanAmax - FanA ) && betaRadian <= pi - 2 * FanA )
% %                                     W = 1 ; 
% %                              elseif ( betaRadian > pi - 2 * FanA && betaRadian <= pi + 2 * FanAmax )        
% %                                     W = sin ( pi * ( pi + 2 * FanAmax - betaRadian ) / ( FanAmax + FanA ) / 4 )^2                                 
% %                              end  
% %                              R ( num ) = W * R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
%                                   R ( num ) = R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% end
% R = reshape ( R , LBeta , LXigama , LP ) ;           % projection data

% % figure , imshow ( squeeze ( R ( : , 331  , : ) ) , [] ) ; 
%% GPU accelerate the preweight & filteration & backprojection
% load R1
z_length = 211;
R= single(R) ;
Display = FDK ( R , Xigamadomain , Pdomain , BetaScanRange , Distance, Size, t_length, s_length, z_length) ;
Display = reshape ( Display , t_length , s_length , z_length ) ;
figure,imshow3Dfull(Display, [0 0.5])
 %% filter    S-L
% load D:\TestCpp\CT\Data\FDK\RFBPimage3\RCT52.mat ;
% load G:\CTcode\Code\ABC.mat
% for i = 1 : 360
%     R( i , : , : ) = Proj ( : , : , i  ) ; 
% end
% 
%  G = zeros ( 2 * LP - 1  , 1) ;                                                      % create filter
% 
% for Pindex = 0 : LP - 1 
%         if Pindex == 0
%             G ( LP ) = 1 / ( pi * PInt )^2 ;    % radian 
%         else                                       
%            G ( LP + Pindex ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * Pindex^2 - 1 ) ;
%            G ( LP - Pindex ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * Pindex^2 - 1 ) ;
%         end
% end
% 
% Rcov = zeros ( LBeta ,  LXigama , LP ) ;
% parfor Betaindex = 1 :  LBeta
%     for Xigamaindex = 1 : LXigama
%     
%         cov = conv ( squeeze ( R ( Betaindex , Xigamaindex , : ) ) , G' ) ;                              % convolution with filter
%         Rcov ( Betaindex , Xigamaindex , : ) = PInt * cov ( LP : 2 * LP - 1 ) ;          
% %         Rcov ( Betaindex , Xigamaindex , : ) = 2 * PInt * cov ( LP : 2 * LP - 1 ) ;                  % half-scan
% 
%     end
% end
% 
% % clear R G ;
% 
% % Hamming = zeros ( length ( Pdomain )  , 1 ) ; 
% % Hamming ( length ( Pdomain ) / 2 - 7 : length ( Pdomain ) / 2 + 7 ) = hamming ( 15 ) ;                   % convolve with hamming window
% % Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% % for i = 1 :  length ( BetaScanRange )
% %         cov = conv ( squeeze ( R ( i , j , : ) ) , Hamming' ) ;                              % convolution with filter
% %         Rcov ( i , j , : ) = cov ( round ( length ( Pdomain ) * 0.5 ) : round ( length ( Pdomain ) * 0.5 ) + length ( Pdomain ) - 1 ) / Hammingsum ;
% % end
% 
% % figure , imshow ( squeeze ( Rcov ( : , 331 , : ) ) , [] ) ; 
% 
% %% reconstruct
% 
% 
% Display = zeros ( size(pic,1) ) ;
% ProjectionSlice = zeros ( size(pic,1) , size(pic,1) , LBeta ) ;
% ProjectionSlice = reshape ( ProjectionSlice , 1 , [] ) ;
% % for t = 1 : t_length
%     
%    t = 257 ;
% %     z = 257 ;
%     parfor num = 1 : z_length * s_length * LBeta
% %                 s = 257 ;
%                z = floor ( ( num - 1 ) / ( LBeta * s_length ) ) + 1 ;
%                
%                s = floor ( ( num - ( z - 1 ) * s_length * LBeta - 1 ) / LBeta ) + 1 ; 
% %                     betas = 1 ;
%                betas = mod ( num - 1 , LBeta ) + 1 ;
%                
%                     beta =  BetaScanRange ( betas ) ;
%                     betaRadian = beta * pi / 180 ;
%                     
%                     source_t = Center_t + Distance * cos ( betaRadian + pi / 2 ) ;                           % define the source
%                     source_s = Center_s + Distance * sin ( betaRadian + pi / 2 ) ;   
%                     source_z = Center_z ; 
%                     
%                     image_t = ( t - 0.5 ) * Resolution2 - Center_t  ;  image_s = ( s - 0.5 ) * Resolution2 - Center_s  ; image_z = ( z - 0.5 ) * Resolution2 - Center_z  ;           % image pixel in ground coordinate
%                     
%                     dect_t = image_t * cos ( betaRadian ) + image_s * sin ( betaRadian ) ;          % in rotate coordinate
%                     dect_s = - image_t * sin ( betaRadian ) + image_s * cos ( betaRadian ) ; 
%                     dect_z = image_z ;     
%                     
%                     LengthRatio = Distance / ( Distance - dect_s ) ; 
%                     Xigama_domain = dect_z * LengthRatio ; 
%                     P_domain = dect_t * LengthRatio ; 
%                     
%                     Xigama_domainN1 =  floor ( ( Xigama_domain + MaxXigama ) / XigamaInt ) + 1 ; 
%                     Xigama_domainN2 = Xigama_domainN1 + 1 ;
%                     P_domainN1 =  floor ( ( P_domain + MaxP ) / PInt ) + 1 ; 
%                     P_domainN2 = P_domainN1 + 1 ;
% 
%                     if  ( P_domain >= -MaxP && P_domain < max ( Pdomain )  && Xigama_domain >= -MaxXigama && Xigama_domain < max ( Xigamadomain ) )
%                             P_domain1 = Pdomain ( P_domainN1 ) ; P_domain2 = Pdomain ( P_domainN2 ) ;
%                             Xigama_domain1 = Xigamadomain ( Xigama_domainN1 ) ; Xigama_domain2 = Xigamadomain ( Xigama_domainN2 ) ;  
%                             
%                             % bilinear interpolation
%                             Xig1 = Xigama_domain - Xigama_domain1 ; Xig2 = Xigama_domain2 - Xigama_domain ;
%                             P1 = P_domain - P_domain1 ; P2 = P_domain2 - P_domain ;
%                             ProjectionSlice ( num ) = ( Xig2 * P2 * Rcov ( betas , Xigama_domainN1 , P_domainN1 ) ...
%                                 + Xig1 * P2 * Rcov ( betas , Xigama_domainN2 , P_domainN1 ) + Xig2 * P1 * Rcov ( betas , Xigama_domainN1 , P_domainN2 ) ...
%                                 + Xig1 * P1 * Rcov ( betas , Xigama_domainN2 , P_domainN2 ) ) / ( PInt * XigamaInt ) * LengthRatio^2 * BetaScanInt ;
%                     end                       
%     end
%      
% % end
% for  z = 1 : z_length 
%     for s = 1 : s_length 
%         Display ( z , s ) = sum ( ProjectionSlice ( 1 + ( ( z - 1 ) * s_length * LBeta + ( s - 1 ) * LBeta ) : LBeta + ( ( z - 1 ) * s_length * LBeta + ( s - 1 ) * LBeta ) ) ) ;
%     end
% end   
% % clear Rcov ProjectionSlice ;
% % figure , imshow ( flipud ( Display ) , [ 0.05 0.35 ] ) ; 
% figure , imshow ( flipud ( Display ) , [ 1 , 1.05 ] ) ; 
% toc ;
% % % figure , plot ( 1 : 513 , squeeze ( Projection ( size ( pic , 1 ) / 2 , : , 50 ) ) , 1 : 513 , squeeze ( pic ( size ( pic , 1 ) / 2 , : , 50 ) )  ) ;
% % % title ( ' grey distrubition ' ) ;
% % % axis ( [ 0 513 1 1.05 ] ) ;
% % % set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;
% % % set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;
% % %% read fig
% % % a = open ( 'D:\TestCpp\CT\Data\FDK\Modified_Para\illed\t257.fig' ) ;
% % % h = get (gca , 'Children') ;
% % % data = get ( h , 'Cdata') ;
% % % figure,imshow ( squeeze ( pic ( 257 , : , : ) )' , [1,1.05] )     s=257
% % 
% % % x=[1:513] ;
% % % y= -8.116e-7 .* ( x -257 ).^2 ;
% % % figure , plot ( 1 : 513  , err ( : , 125 ) , 1 : 513 , y );axis ( [ 0 513 -0.05 0 ] ) ;

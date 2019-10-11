%2016/11/28 ZXZ
% equal spaced
% reference : Algorithms for Reconstruction with Nondiffracting Sources
%              |                     coordinate system in matlab different from the one in the theory
%              |
%              |
% -----------|------------> j
%              |
%              |
%             \/  i
%  **** RScan between source and center point play a great role in suppressing unintersted area *****
%  
%  % we define the actual length of the pic as ( x = 60 mm y = 60 mm )
%   note that resolution is different from length

clear all;
Size = [ 60 , 60 ] ;                                  % actual range
pic = single(phantom( 256 )) ;           % original picture : number means the number of pixel in the range 

BetaScanInt = deg2rad(0.5) ;             % scanning internal              
MaxBeta = deg2rad(360) ; 
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project

MaxP = RPic * ( 1 + 0.1 )  ;                                           
PInt = 0.1 ;                      %   interval of S ( interval on the detect plain )
Pdomain = - MaxP : PInt : MaxP ;                          % detective range
LP = length ( Pdomain ) ;

Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  

Ratio = 4 ;                                                           % shoul d be smaller than 8
RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

R = zeros ( LBeta ,  LP ) ;   % create space to store fan projection
picvector = Img2vec_Mat2Cpp2D( pic ) ;
%%  GPU projection
R = ProjectionFan_2D ( picvector, height, width, Size, BetaScanRange', Pdomain', RScan ) ;
% R = reshape( R , LP , LBeta )' ;
% figure,imshow(R, [])
%% system matrix
% SysMatrix = GenSysMatFan ( height, width, Size, BetaScanRange, Pdomain, RScan, Center_x , Center_y) ;
% R = SysMatrix * double(picvector) ;
R = reshape( R , LP , LBeta )' ;
% figure,imshow(R, [])

%%  fan derived from inner radon projecion
% for test step one 
% thetaInt = 1 ; 
% theta = 1 : thetaInt : 360 ; 
% [ Radon , xp ] = radon ( pic , theta ) ;
% for i = 1 :  length ( BetaScanRange )
%      beta =  BetaScanRange ( i ) ; 
%      for j = 1 :  length ( Pdomain )
%          gamaRadian = atan ( Pdomain ( j ) / RScan ) ;
%          gama = gamaRadian * 180 / pi ;
%          t_axis = RScan * sin ( gamaRadian ) + 0.5 * size ( Radon , 1 ) + 1 ; 
%          theta_axis = round ( beta + gama ) ; 
%          if ( theta_axis > 360 )
%                 theta_axis = theta_axis - 360 ; 
%          elseif ( theta_axis <= 0 )
%                 theta_axis = theta_axis + 360 ;  
%          end
%          R ( i , j ) = Radon ( round ( t_axis ) , theta_axis ) ;
%          R ( i , j ) = R ( i , j ) * cos ( gamaRadian ) ;
%      end
% end
% figure , imagesc( R ) ; 

%% direct fan

%step one 
% RScanInt = 1 ;      %define the interval on the rayline
% for i = 1 :  length ( BetaScanRange ) 
%         beta =  BetaScanRange ( i ) ;
%         betaRadian = beta * pi / 180 ;                    % angle to radian
%         source_i = Center_i + RScan * cos ( betaRadian + pi / 2 ) ;  source_j = Center_i + RScan * sin ( betaRadian + pi / 2 ) ;  %define the source in matlab coordinate
%         for j = 1 : length ( Pdomain )
%                     gamaRadian = atan ( Pdomain ( j ) / RScan ) ;
% 
%                     Smax = RScan + 2 * RPic ; 
%                     DetectPoint_iend = Center_i + Smax * sin ( gamaRadian + betaRadian ) - RScan * sin ( betaRadian ) ;   % define end detect point in matlab coordinate
%                     DetectPoint_jend = Center_i - Smax * cos ( gamaRadian + betaRadian ) + RScan * cos ( betaRadian ) ;
%                     if ( floor ( DetectPoint_iend ) == floor ( source_i ) || ceil ( DetectPoint_iend ) == ceil ( source_i ) )                     % define projection line direction and range of i and j
%                          i_range = [ ] ;
%                     else
%                         if ( DetectPoint_iend > source_i )
%                              i_signal = 1 ;   
%                              if ( source_i < 0 )
%                                      istart = 0 ;
%                              else
%                                      istart = ceil ( source_i / Resolution ) ;
%                              end
%                              if ( DetectPoint_iend > Size ( 1 ) )
%                                      iend = height ;
%                              else
%                                      iend = floor ( DetectPoint_iend / Resolution ) ;
%                              end                         
%                         elseif ( DetectPoint_iend < source_i )
%                             i_signal = -1 ; 
%                             if ( DetectPoint_iend < 0 )
%                                      istart = 0 ;
%                              else
%                                      istart = ceil ( DetectPoint_iend / Resolution ) ;
%                              end
%                              if ( source_i > Size ( 1 ) )
%                                      iend = height ;
%                              else
%                                      iend = floor ( source_i / Resolution ) ;
%                              end
%                         end
%                         i_range = istart : iend ;
%                     end
%                     
%                     if ( floor ( DetectPoint_jend ) == floor ( source_j ) || ceil ( DetectPoint_jend ) == ceil ( source_j ) )
%                          j_range = [ ] ;
%                     else
%                         if ( DetectPoint_jend > source_j )
%                             j_signal = 1 ; 
%                             if ( source_j < 0 )
%                                      jstart = 0 ;
%                              else
%                                      jstart = ceil ( source_j / Resolution ) ;
%                              end
%                              if ( DetectPoint_jend > Size ( 2 ) )
%                                      jend = width ;
%                              else
%                                      jend = floor ( DetectPoint_jend / Resolution ) ;
%                              end
%                         elseif ( DetectPoint_jend < source_j )
%                             j_signal = -1 ; 
%                             if ( DetectPoint_jend < 0 )
%                                      jstart = 0 ;
%                              else
%                                      jstart = ceil ( DetectPoint_jend / Resolution ) ;
%                              end
%                              if ( source_j > Size ( 2 ) )
%                                      jend = width ;
%                              else
%                                      jend = floor ( source_j / Resolution ) ;
%                              end
%                         end
%                         j_range = jstart : jend ;
%                     end
%                     
%                     Li = length ( i_range ) ; Lj = length ( j_range ) ;
%                     ProjectionLine = zeros ( 3 , Li  + Lj ) ;
%                     ProjectionLine ( 1 , 1 : Li  ) = abs ( ( i_range * Resolution - source_i ) / ( DetectPoint_iend - source_i ) ) * Smax ;
%                     ProjectionLine ( 2 , 1 : Li  ) = i_range ;
%                     ProjectionLine ( 3 , 1 : Li  ) = 1 ;                              % 1 represents i
%                     ProjectionLine ( 1 , ( Li  + 1 ) : ( Li  + Lj ) ) = abs ( ( j_range * Resolution - source_j ) / ( DetectPoint_jend - source_j ) ) * Smax ;
%                     ProjectionLine ( 2 , ( Li  + 1 ) : ( Li  + Lj ) ) = j_range ;
%                     ProjectionLine ( 3 , ( Li  + 1 ) : ( Li  + Lj ) ) = 2 ;                              % 2 represents j
%                     ProjectionLine = sortrows ( ProjectionLine' )' ;
%                     DetectPoint_i = ceil ( source_i / Resolution ) ; DetectPoint_j = ceil ( source_j / Resolution ) ;                   % start point of the line        
% 
%                     for n = 1 : length ( ProjectionLine )
% 
%                         if ( ProjectionLine ( 3 , n ) == 1 && Li ~= 0 )                                                              % define the pixel 
%                                 DetectPoint_i = ProjectionLine ( 2 , n ) + 0.5 + i_signal * 0.5 ; 
%                         elseif ( ProjectionLine ( 3 , n ) == 2 && Lj ~= 0 )
%                                 DetectPoint_j = ProjectionLine ( 2 , n ) + 0.5 + j_signal * 0.5 ; 
%                         end
% 
%                         if ( DetectPoint_i > 0 && DetectPoint_i <= height && DetectPoint_j > 0 && DetectPoint_j <= width ) 
%                                  f = pic ( DetectPoint_i , DetectPoint_j ) ;
%                                  LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                                  R ( i , j ) = R ( i , j ) + f * LengthUnit ;    
%                          end 
%                    end 
% 
%                    R ( i , j ) = R ( i , j ) * cos ( gamaRadian ) ;
%         end
% end
% figure , imshow ( R , [] ) ; 
%% projection in formula  P.54
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
%      
% proportion = max ( Size ) / 2 ;                            % actual scale
% R = reshape ( R , 1 , [] ) ; 
% 
% for num = 1 : LBeta * LP 
%     
%         Betaindex = mod ( num - 1 ,  LBeta ) + 1 ;
%         Sindex = floor ( ( num - 1 ) / LBeta ) + 1 ;
%         
%         beta =  BetaScanRange ( Betaindex ) ;
%         betaRadian = beta * pi / 180 ;                    % angle to radian
%         %define the source in matlab coordinate
%         source_x = Center_x + RScan * cos ( betaRadian + pi / 2 ) ;  source_y = Center_y + RScan * sin ( betaRadian + pi / 2 ) ;         % in matlab coordinate
%                    
%                     gamaRadian = atan ( Pdomain ( Sindex ) / RScan ) ;
%                     t = RScan * sin ( gamaRadian ) ;
%                     thetaRadian = gamaRadian + betaRadian ;
%                     for n = 1 : 10
%                         x0 = shep ( n , 4 ) * proportion ; y0 = shep ( n , 5 ) * proportion ;              % oval information in ground coordinate
%                         s = sqrt ( x0^2 + y0^2 ) ;
%                         ggama = atan2 ( y0 , x0 ) ;
%                         alpaha = shep ( n , 6 ) * pi / 180 ;
%                         A = shep ( n , 2 ) * proportion ; B = shep ( n , 3 ) * proportion  ;    
%                         Rou = shep ( n , 1 ) ;                   
%                         ts = t - s * cos ( ggama - thetaRadian ) ;
%                         thetas = thetaRadian - alpaha ; 
%                         a2 = ( A * cos ( thetas ) )^2 + ( B * sin ( thetas ) )^2 ;
%                         if ( abs ( ts ) <= sqrt ( a2 ) )                                             % decide whether cross the oval
%                                 R ( num ) = R ( num ) + 2 * Rou * A * B * sqrt ( a2 - ts^2 ) / a2 ;
%                         end                             
%                     end   
%                     R ( num ) = R ( num ) * cos ( gamaRadian ) ;
% end
% R = reshape ( R , LBeta , LP ) ;
% figure , imshow ( R , [] ) ; 
%  for i = 1 : LP
%     N = i - ( LP + 1 ) / 2 ;
%         if N == 0
%             G ( i ) = 1 / ( 8 * PInt ^2 ) ;    % radian 
%         elseif rem ( N , 2 ) == 0          %even
%             G ( i ) = 0 ;
%         else                                         %odd
%            G ( i ) = - 0.5 / ( pi * Pdomain ( i ) )  ^ 2 ;
%         end
% end
%% filter convolution based R-L
%%step two
%
% Rcov = zeros ( length ( BetaScanRange ) ,  length ( Pdomain ) ) ;
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = PInt * cov ( round ( length ( Pdomain ) * 0.5 ) : round ( length ( Pdomain ) * 0.5 ) + length ( Pdomain ) - 1 ) ;
%        % continuous domain to discret domain
%         % Rcov ( i , : ) = cov ( 1 :  length ( GamaDetectRange )  )  ;
%         %Rcov ( i , : ) = cov (  length ( GamaDetectRange ) : 2 * length ( GamaDetectRange ) - 1 )  ;
% end
% G = zeros ( 2 * length ( Pdomain ) - 1  , 1) ;                                            % create filter
% LP = length ( Pdomain ) ; 
% for i = 0 : LP - 1 
%         if i == 0
%             G ( LP ) = 1 / ( 8 * PInt ^2 ) ;    % radian 
%         elseif rem ( i , 2 ) == 0          % even
%             G ( LP + i ) = 0 ;
%             G ( LP - i ) = 0 ;
%         else                                         % odd
%            G ( LP + i ) = - 0.5 / ( pi * i * PInt ) ^ 2 ;
%            G ( LP - i ) = - 0.5 / ( pi * i * PInt  )  ^ 2 ;
%         end
% end
% LBeta = length ( BetaScanRange ) ; 
% Rcov = zeros ( LBeta ,  LP ) ;
% for i = 1 :  LBeta
%         cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = PInt * cov ( LP : 2 * LP - 1 ) ;
% end
% 
% Hamming = zeros ( 2 * LP - 1 , 1 ) ; 
% Hamming ( LP - 5 : LP + 5 ) = hamming ( 11 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( Rcov ( i , : ) , Hamming' ) ;                              % convolution with filter
%         Rcov ( i , : ) = cov ( LP : 2 * LP - 1 ) / Hammingsum ;
% end

% figure , imshow ( Rcov ) ; 

%% filter convolution based S-L
%%step two
tic 
G = zeros ( 2 * LP - 1  , 1) ;                                            % create filter

for i = 0 : LP - 1 
        if i == 0
            G ( LP ) = 1 / ( pi * PInt )^2 ;    % radian 
        else                                       
           G ( LP + i ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * i^2 - 1 ) ;
           G ( LP - i ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * i^2 - 1 ) ;
        end
end

Rcov = zeros ( LBeta ,  LP ) ;
parfor i = 1 :  LBeta
        cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
        Rcov ( i , : ) = PInt * cov ( LP : 2 * LP - 1 ) ;
end
clear R G ;
% figure , imagesc( Rcov ) ; 
%% filt in frequency domain
% LP = length ( Pdomain ) ; 
% G = zeros ( 1 , LP ) ;      
% LGh = round ( ( LP  - 1 ) / 2 ) ;
% for i = 0 : LGh - 1
%         if i == 0
%             G ( LGh ) = 0 ;    % radian 
%         else                                       
%            G ( LGh + i ) = i / ( 2 * PInt * ( LGh -1 ) ) ;
%            G ( LGh - i ) = i / ( 2 * PInt * ( LGh -1 ) ) ;
%         end  
% end
% LBeta = length ( BetaScanRange ) ; 
% Rcov = zeros ( LBeta ,  LP ) ;
% for i = 1 :  LBeta
%         Rfft = fftshift ( fft ( R ( i , : ) ) ) ;                                             
%         Rm = Rfft .* G ;                         % convolution with filter
%         Rcov ( i , : ) = real ( ifft ( ifftshift ( Rm ) ) ) ;
% end

%%  reconstruct
%step three
% 
ProjectionShift = zeros ( LBeta , height , width ) ;                           % temporarily store back-porjection
ProjectionShift = reshape ( ProjectionShift , 1 , [] ) ;
Projection = zeros ( height , width ) ;                                    %  reconsruct picture of 513*513 resolution
Resolution2 =  max ( Size ) / height ;                                 % define the resolution of reconstruction

parfor num = 1 :  LBeta * height * width 
                                                 
                    i = floor ( ( num - 1 ) / ( height * LBeta ) ) + 1 ;                % piexl number
                    j = floor ( ( num - ( i - 1 ) * width * LBeta - 1 )  / LBeta ) + 1 ; 
                    beta_index = mod ( num - 1 , LBeta ) + 1 ;                  
                    
                    Pixel_x = ( i - 0.5 ) * Resolution2 - Center_x ;        % pixel in ground coordinate
                    Pixel_y = ( j - 0.5 ) * Resolution2 - Center_y;   
                    image_pixel = [ Pixel_x ; Pixel_y ] ;           
                    
                    beta =  BetaScanRange ( beta_index ) ;
%                     source_x = RScan * cos ( betaRadian + pi / 2 ) ;  source_y = RScan * sin ( betaRadian + pi / 2 ) ;  %define the source in ground coordinate
%                     source_pixel = [ source_x , source_y ] ;
%                
                    H = zeros ( 2 ) ;                                                  % rotate matrix
                    H ( 1 , 1 ) = cos ( beta ) ;
                    H ( 1 , 2 ) = sin ( beta ) ;
                    H ( 2 , 1 ) = -sin ( beta ) ;
                    H ( 2, 2 ) = cos ( beta ) ;                   
                    Image_m = H * image_pixel ;
                    
                    U = ( RScan - Image_m ( 2 ) ) / RScan ;
                    Image_shiftpixel = [ 1 0 ; 0 -1 ] * ( Image_m - [ 0 ; RScan ] ) ;
                    gama_shift = atan ( Image_shiftpixel ( 1 ) / Image_shiftpixel ( 2 ) ) ;
                    
                    S_domain = RScan * Image_shiftpixel ( 1 ) / Image_shiftpixel ( 2 ) ;                        
                    S_domainN1 = floor ( ( S_domain + MaxP ) / PInt ) + 1 ;
                    S_domainN2 = S_domainN1 + 1 ;
                    
                    if  ( S_domain >= min ( Pdomain ) && S_domain < max ( Pdomain ) )
                            S1 =  S_domain - Pdomain ( S_domainN1 ) ;
                            S2 =  Pdomain ( S_domainN2 ) - S_domain ;
                            ProjectionShift ( num ) = ( S2 * Rcov ( beta_index , S_domainN1 ) + ...
                                S1 * Rcov ( beta_index , S_domainN2 ) ) / ( U^2 * PInt ) * BetaScanInt ;
                    end                        
end
clear Rocv ;

for i = 1 : width 
    for j = 1 : height 
            Projection ( i , j ) = sum ( ProjectionShift ( 1 + ( i - 1 ) * LBeta * height + ( j - 1 ) * LBeta : LBeta + ( i - 1 ) * LBeta * width + ( j - 1 ) * LBeta ) ) ;
    end
end

toc
%% Display 

ProjectionDisplay = flipud ( Projection' ) ;
figure , imshow( ProjectionDisplay , [ 0 , 1 ] ) ; 
title ( ' Reconstructed image ' ) ;
pic = phantom( height , width ) ;
mid_index = round(height/2) ;
figure , plot ( 1 : size ( pic , 1 ) , ProjectionDisplay ( mid_index , : ) , 1 : height , pic ( mid_index , : ) ) ;
title ( ' grey distrubition ' ) ;
axis ( [ 0 height 0 1 ] ) ;
%%   frequency display
% Projectionfft = fft2 ( Projection ) ;
% Projectionffts = fftshift ( Projectionfft ) ;
% Fm = abs ( Projectionffts ) ;  
% figure, imshow ( log ( 1  + Fm ) , [ ] )

%% evaluation function
% aver=sum ( sum ( pic ) ) / ( size ( pic , 1 ) * size ( pic , 2 ) ) ;
% pd=double ( pic ) ;
% d= ( sum ( sum ( ( pd - Projection ).^2 ) ) / sum ( sum ( ( pd - aver ).^2 ) ) ) ^0.5 ;

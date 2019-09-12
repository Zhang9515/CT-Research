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
tic 
clear all;
Size = [ 60 , 60 ] ;                                  % actual range
pic = StandardPhantom( 513 ) ;           % original picture : number means the number of pixel in the range 

BetaScanInt = 1 ;             % scanning internal              
MaxBeta = 360 ; 
BetaScanRange = BetaScanInt + 0.1 : BetaScanInt : MaxBeta + 0.1  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

[ HeightN , WidthN ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / max ( size ( pic ) ) ;   % define the resolution of pic
RPic = max ( Size ) * sqrt ( 2 ) / 2 ;                     % radius of project

MaxS = RPic * ( 1 + 0.1 )  ;                                           
SInt = 0.1 ;                      %   interval of S ( interval on the detect plain )
Sdomain = - MaxS : SInt : MaxS ;                          % detective range
LS = length ( Sdomain ) ;

Center_x = Size ( 1 ) / 2 ;  Center_y = Size ( 2 ) / 2 ;      % make the center point overlay the center pixel  

Ratio = 4 ;                                                           % shoul d be smaller than 8
RScan = RPic * Ratio ;                                        % distance between source and center point ( radius of trajectory ) 

R = zeros ( length ( BetaScanRange ) ,  length ( Sdomain ) ) ;   % create space to store fan projection
%%  GPU projection
picvector = reshape (pic, HeightN * WidthN, 1);
R = ProjectionFan_2D (picvector, HeightN, WidthN, Size, BetaScanRange', Sdomain', RScan);
R = reshape( R , LS , LBeta ) ;
figure,imshow(R, [])

%%  fan derived from inner radon projecion
% for test step one 
% thetaInt = 1 ; 
% theta = 1 : thetaInt : 360 ; 
% [ Radon , xp ] = radon ( pic , theta ) ;
% for i = 1 :  length ( BetaScanRange )
%      beta =  BetaScanRange ( i ) ; 
%      for j = 1 :  length ( Sdomain )
%          gamaRadian = atan ( Sdomain ( j ) / RScan ) ;
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
%         for j = 1 : length ( Sdomain )
%                     gamaRadian = atan ( Sdomain ( j ) / RScan ) ;
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
%                                      iend = HeightN ;
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
%                                      iend = HeightN ;
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
%                                      jend = WidthN ;
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
%                                      jend = WidthN ;
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
%                         if ( DetectPoint_i > 0 && DetectPoint_i <= HeightN && DetectPoint_j > 0 && DetectPoint_j <= WidthN ) 
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
     
proportion = max ( Size ) / 2 ;                            % actual scale
R = reshape ( R , 1 , [] ) ; 

for num = 1 : LBeta * LS 
    
        Betaindex = mod ( num - 1 ,  LBeta ) + 1 ;
        Sindex = floor ( ( num - 1 ) / LBeta ) + 1 ;
        
        beta =  BetaScanRange ( Betaindex ) ;
        betaRadian = beta * pi / 180 ;                    % angle to radian
        %define the source in matlab coordinate
        source_x = Center_x + RScan * cos ( betaRadian + pi / 2 ) ;  source_y = Center_y + RScan * sin ( betaRadian + pi / 2 ) ;         % in matlab coordinate
                   
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
                    R ( num ) = R ( num ) * cos ( gamaRadian ) ;
end
R = reshape ( R , LBeta , LS ) ;
% figure , imshow ( R , [] ) ; 
 for i = 1 : length ( Sdomain )
    N = i - ( length ( Sdomain ) + 1 ) / 2 ;
        if N == 0
            G ( i ) = 1 / ( 8 * SInt ^2 ) ;    % radian 
        elseif rem ( N , 2 ) == 0          %even
            G ( i ) = 0 ;
        else                                         %odd
           G ( i ) = - 0.5 / ( pi * Sdomain ( i ) )  ^ 2 ;
        end
end
%% filter convolution based R-L
%%step two
%
% Rcov = zeros ( length ( BetaScanRange ) ,  length ( Sdomain ) ) ;
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = SInt * cov ( round ( length ( Sdomain ) * 0.5 ) : round ( length ( Sdomain ) * 0.5 ) + length ( Sdomain ) - 1 ) ;
%        % continuous domain to discret domain
%         % Rcov ( i , : ) = cov ( 1 :  length ( GamaDetectRange )  )  ;
%         %Rcov ( i , : ) = cov (  length ( GamaDetectRange ) : 2 * length ( GamaDetectRange ) - 1 )  ;
% end
% G = zeros ( 2 * length ( Sdomain ) - 1  , 1) ;                                            % create filter
% Ls = length ( Sdomain ) ; 
% for i = 0 : Ls - 1 
%         if i == 0
%             G ( Ls ) = 1 / ( 8 * SInt ^2 ) ;    % radian 
%         elseif rem ( i , 2 ) == 0          % even
%             G ( Ls + i ) = 0 ;
%             G ( Ls - i ) = 0 ;
%         else                                         % odd
%            G ( Ls + i ) = - 0.5 / ( pi * i * SInt ) ^ 2 ;
%            G ( Ls - i ) = - 0.5 / ( pi * i * SInt  )  ^ 2 ;
%         end
% end
% LBeta = length ( BetaScanRange ) ; 
% Rcov = zeros ( LBeta ,  Ls ) ;
% for i = 1 :  LBeta
%         cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = SInt * cov ( Ls : 2 * Ls - 1 ) ;
% end
% 
% Hamming = zeros ( 2 * Ls - 1 , 1 ) ; 
% Hamming ( Ls - 5 : Ls + 5 ) = hamming ( 11 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( Rcov ( i , : ) , Hamming' ) ;                              % convolution with filter
%         Rcov ( i , : ) = cov ( Ls : 2 * Ls - 1 ) / Hammingsum ;
% end

% figure , imshow ( Rcov ) ; 

%% filter convolution based S-L
%%step two

G = zeros ( 2 * LS - 1  , 1) ;                                            % create filter

for i = 0 : LS - 1 
        if i == 0
            G ( LS ) = 1 / ( pi * SInt )^2 ;    % radian 
        else                                       
           G ( LS + i ) = - 1 / ( pi * SInt ) ^ 2 / ( 4 * i^2 - 1 ) ;
           G ( LS - i ) = - 1 / ( pi * SInt ) ^ 2 / ( 4 * i^2 - 1 ) ;
        end
end

Rcov = zeros ( LBeta ,  LS ) ;
parfor i = 1 :  LBeta
        cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
        Rcov ( i , : ) = SInt * cov ( LS : 2 * LS - 1 ) ;
end
clear R G ;
% figure , imagesc( Rcov ) ; 
%% filt in frequency domain
% Ls = length ( Sdomain ) ; 
% G = zeros ( 1 , Ls ) ;      
% LGh = round ( ( Ls  - 1 ) / 2 ) ;
% for i = 0 : LGh - 1
%         if i == 0
%             G ( LGh ) = 0 ;    % radian 
%         else                                       
%            G ( LGh + i ) = i / ( 2 * SInt * ( LGh -1 ) ) ;
%            G ( LGh - i ) = i / ( 2 * SInt * ( LGh -1 ) ) ;
%         end  
% end
% LBeta = length ( BetaScanRange ) ; 
% Rcov = zeros ( LBeta ,  Ls ) ;
% for i = 1 :  LBeta
%         Rfft = fftshift ( fft ( R ( i , : ) ) ) ;                                             
%         Rm = Rfft .* G ;                         % convolution with filter
%         Rcov ( i , : ) = real ( ifft ( ifftshift ( Rm ) ) ) ;
% end

%%  reconstruct
%step three
% 
ProjectionShift = zeros ( LBeta , 513 , 513 ) ;                           % temporarily store back-porjection
ProjectionShift = reshape ( ProjectionShift , 1 , [] ) ;
Projection = zeros ( 513 , 513 ) ;                                    %  reconsruct picture of 513*513 resolution
Resolution2 =  max ( Size ) / 513 ;                                 % define the resolution of reconstruction

parfor num = 1 :  LBeta * 513 * 513 
                                                 
                    i = floor ( ( num - 1 ) / ( 513 * LBeta ) ) + 1 ;                % piexl number
                    j = floor ( ( num - ( i - 1 ) * 513 * LBeta - 1 )  / LBeta ) + 1 ; 
                    betas = mod ( num - 1 , LBeta ) + 1 ;                  
                    
                    Pixel_x = ( i - 0.5 ) * Resolution2 - Center_x ;        % pixel in ground coordinate
                    Pixel_y = ( j - 0.5 ) * Resolution2 - Center_y;   
                    image_pixel = [ Pixel_x ; Pixel_y ] ;           
                    
                    beta =  BetaScanRange ( betas ) ;
                    betaRadian = beta * pi / 180 ;
%                     source_x = RScan * cos ( betaRadian + pi / 2 ) ;  source_y = RScan * sin ( betaRadian + pi / 2 ) ;  %define the source in ground coordinate
%                     source_pixel = [ source_x , source_y ] ;
%                
                    H = zeros ( 2 ) ;                                                  % rotate matrix
                    H ( 1 , 1 ) = cos ( betaRadian ) ;
                    H ( 1 , 2 ) = sin ( betaRadian ) ;
                    H ( 2 , 1 ) = -sin ( betaRadian ) ;
                    H ( 2, 2 ) = cos ( betaRadian ) ;                   
                    Image_m = H * image_pixel ;
                    
                    U = ( RScan - Image_m ( 2 ) ) / RScan ;
                    Image_shiftpixel = [ 1 0 ; 0 -1 ] * ( Image_m - [ 0 ; RScan ] ) ;
                    gama_shift = atan ( Image_shiftpixel ( 1 ) / Image_shiftpixel ( 2 ) ) ;
                    
                    S_domain = RScan * Image_shiftpixel ( 1 ) / Image_shiftpixel ( 2 ) ;                        
                    S_domainN1 = floor ( ( S_domain + MaxS ) / SInt ) + 1 ;
                    S_domainN2 = S_domainN1 + 1 ;
                    
                    if  ( S_domain >= min ( Sdomain ) && S_domain < max ( Sdomain ) )
                            S1 =  S_domain - Sdomain ( S_domainN1 ) ;
                            S2 =  Sdomain ( S_domainN2 ) - S_domain ;
                            ProjectionShift ( num ) = ( S2 * Rcov ( betas , S_domainN1 ) + ...
                                S1 * Rcov ( betas , S_domainN2 ) ) / ( U^2 * SInt ) * BetaScanInt * pi / 180 ;
                    end                        
end
clear Rocv ;

for i = 1 : 513 
    for j = 1 : 513 
            Projection ( i , j ) = sum ( ProjectionShift ( 1 + ( i - 1 ) * LBeta * 513 + ( j - 1 ) * LBeta : LBeta + ( i - 1 ) * LBeta * 513 + ( j - 1 ) * LBeta ) ) ;
    end
end


%% Display 

ProjectionDisplay = flipud ( Projection' ) ;
figure , imshow( ProjectionDisplay , [ 1 , 1.05 ] ) ; 
title ( ' Reconstructed image ' ) ;
pic = StandardPhantom( 513 ) ;
figure , plot ( 1 : size ( pic , 1 ) , ProjectionDisplay ( 229 , : ) , 1 : 513 , pic ( 229 , : ) ) ;
title ( ' grey distrubition ' ) ;
axis ( [ 0 513 1 1.05 ] ) ;
%%   frequency display
% Projectionfft = fft2 ( Projection ) ;
% Projectionffts = fftshift ( Projectionfft ) ;
% Fm = abs ( Projectionffts ) ;  
% figure, imshow ( log ( 1  + Fm ) , [ ] )

%% evaluation function
% aver=sum ( sum ( pic ) ) / ( size ( pic , 1 ) * size ( pic , 2 ) ) ;
% pd=double ( pic ) ;
% d= ( sum ( sum ( ( pd - Projection ).^2 ) ) / sum ( sum ( ( pd - aver ).^2 ) ) ) ^0.5 ;
Time = toc ;
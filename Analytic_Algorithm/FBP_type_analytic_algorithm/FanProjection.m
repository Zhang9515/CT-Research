%2016/11/16 ZXZ
% equiangular ray
% reference : Algorithms for Reconstruction with Nondiffracting Sources
%              |                     coordinate system in matlab
%              |
% -----------|------------> j
%              |
%              |
%             \/  i
tic
clear all;
BetaScanInt = 240 / 800 * pi / 180 ;             % scanning internal
GamaDetectInt = 9.7527e-4 ;              % detective internal                 
% pic = phantom( 512 ) ;           % original picture
%pic=dicomread('1.3.6.1.4.1.9328.50.6.579.dcm');
%figure , imagesc( pic ) ;             % show in fake colour
MaxBeta = 240 * pi / 180 ; 
StartAngle = 1.5466 * pi / 180 ; 
BetaScanRange = StartAngle : BetaScanInt : ( MaxBeta + StartAngle - BetaScanInt ) ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

% [ height , width ] = size ( pic ) ;              % store the size of picture
height  = 512 ; width = 512 ;
Center_i = height / 2 ;  Center_j = width / 2 ;      % define the center 
% Distance = max ( size ( pic ) ) ;                        % define the distance between center and source
Distance = 570 ;
% MaxGama = round ( atan ( 0.5 * sqrt ( 2 ) * height / Distance ) * 180 / pi ) ; 

% GamaDetectRange = -MaxGama : GamaDetectInt : MaxGama ;    % detective range
LGama =  936 ;                                                                                      % Gama ( radian ) 
fid = fopen( 'G:\CTcode\Data\UIData\1\detectorangle.dat' , 'rb' ) ;         
[ GamaDetectRange ,  ~ ] = fread ( fid , LGama , 'float' ) ;
% GamaDetectRange = GamaDetectRange - 0.004876 ; % modification
GamaDetectRange = -GamaDetectRange ;
fclose ( fid ) ;
MaxGama = abs ( GamaDetectRange ( end ) - GamaDetectRange ( 1 ) ) ; %%%%%

% R = zeros ( LBeta ,  LGama ) ;   % create space to store fan projection
fid = fopen ( 'G:\CTcode\Data\UIData\1\data.dat' , 'rb' ) ;                                   % load projection data
[ Rs , count ] = fread ( fid , LGama * 80 * 800 , 'float' ) ;
Rs = reshape ( Rs , LGama , 80 , 800 ) ;
fclose ( fid ) ;
R = squeeze ( Rs ( : , 40 , : ) ) ; 
% R = flipud ( R )' ;       % inverse the Gama domain
R = R' ;

G = zeros ( LGama  , 1 ) ;                                            % create filter

Projection = zeros ( height , width ) ;      %reconstruct the pic

%%  fan derived from inner radon projecion
% for test step one 
% thetaInt = 1 ; 
% theta = 1 : thetaInt : 360 ; 
% [ Radon , xp ] = radon ( pic , theta ) ;
% for i = 1 :  length ( BetaScanRange )
%      beta =  BetaScanRange ( i ) ; 
%      for j = 1 : length ( GamaDetectRange )
%          gama =  GamaDetectRange ( j ) ; 
%          gamaRadian = gama * pi / 180 ;
%          t_axis = Distance * sin ( gamaRadian ) + 0.5 * size ( Radon , 1 ) + 1 ; 
%          theta_axis = round ( beta + gama ) ; 
%          if ( theta_axis > 360 )
%                 theta_axis = theta_axis - 360 ; 
%          elseif ( theta_axis <= 0 )
%                 theta_axis = theta_axis + 360 ;  
%          end
%          R ( i , j ) = Radon ( round ( t_axis ) , theta_axis ) ;
%          R ( i , j ) = R ( i , j ) * Distance * cos ( gamaRadian ) ;
%      end
% end
% figure , imagesc( R ) ; 

%% fan derived from my radon projecion
% thetaInt = 1 ; 
% theta = 1 : thetaInt : 360 ;
% tmax = round ( 0.5 * sqrt ( 2 ) * size ( pic , 1 ) ) ;
% Radon = zeros ( 2 * tmax + 1 ,  length ( theta ) ) ;   % create space to store fan projection
% distanceint = 1 ; 
% 
% for theta_axis = theta
%     t_axis = 0 ;
%     while t_axis <= tmax 
%         dis = 0 ;
%         dect_i = Center_i - t_axis * sin ( theta_axis * pi / 180 ) ; dect_j = Center_j + t_axis * cos ( theta_axis * pi / 180 ) ; 
%         pixel_i = dect_i ; pixel_j = dect_j ; pixel_ii = dect_i ; pixel_jj = dect_j ;
%         while ( ( pixel_i > 0 && pixel_i < height && pixel_j > 0 && pixel_j < width ) || ( pixel_ii > 0 && pixel_ii < height && pixel_jj > 0 && pixel_jj < width ) )
%             Radon ( t_axis + tmax + 1 , theta_axis ) = Radon ( t_axis + tmax + 1 , theta_axis ) + distanceint * pic ( floor ( pixel_i ) + 1 , floor ( pixel_j ) + 1 ) ;           
%             if ( dis ~= 0 )
%                 Radon ( t_axis + tmax + 1 , theta_axis ) = Radon ( t_axis + tmax + 1 , theta_axis ) + distanceint * pic ( floor ( pixel_ii ) +1 , floor ( pixel_jj ) + 1 ) ;
%             end 
%             dis = dis + distanceint ;  
%             pixel_i = dect_i - dis * cos ( theta_axis * pi / 180 ) ; pixel_j = dect_j - dis * sin ( theta_axis * pi / 180 ) ;     
%             pixel_ii = dect_i + dis * cos ( ( theta_axis + 180 ) * pi / 180 ) ; pixel_jj = dect_j + dis * sin ( ( theta_axis + 180 ) * pi / 180 ) ;
%         end
%         t_axis = t_axis + 1 ;
%     end 
%     
%     t_axis = -1 ;
%     while t_axis >= -tmax 
%         dis = 0 ;
%         dect_i = Center_i - t_axis * sin ( theta_axis * pi / 180 ) ; dect_j = Center_j + t_axis * cos ( theta_axis * pi / 180 ); 
%         pixel_i = dect_i ; pixel_j = dect_j ; pixel_ii = dect_i ; pixel_jj = dect_j ;
%          while ( ( pixel_i > 0 && pixel_i < height && pixel_j > 0 && pixel_j < width ) || ( pixel_ii > 0  && pixel_ii < height && pixel_jj > 0 && pixel_jj < width ) )
%             Radon ( t_axis + tmax + 1 , theta_axis ) = Radon ( t_axis + tmax + 1 , theta_axis ) + distanceint * pic ( floor ( pixel_i ) + 1 , floor ( pixel_j ) + 1 ) ;
%             if ( dis ~= 0 )
%                 Radon ( t_axis + tmax + 1 , theta_axis ) = Radon ( t_axis + tmax + 1 , theta_axis ) + distanceint * pic ( floor ( pixel_ii ) + 1 , floor ( pixel_jj ) + 1 ) ;
%             end 
%             dis = dis + distanceint ;     
%             pixel_i = dect_i - dis * cos ( theta_axis * pi / 180 ) ; pixel_j = dect_j - dis * sin ( theta_axis * pi / 180 ) ;      
%             pixel_ii = dect_i + dis * cos ( ( theta_axis + 180 ) * pi / 180 ) ; pixel_jj = dect_j + dis * sin ( ( theta_axis + 180 ) * pi / 180 ) ;
%         end
%         t_axis = t_axis - 1 ;
%     end 
% end
% figure,imagesc(Radon); 
% 
% for i = 1 :  length ( BetaScanRange )
%      beta =  BetaScanRange ( i ) ; 
%      for j = 1 : length ( GamaDetectRange )
%          gama =  GamaDetectRange ( j ) ; 
%          gamaRadian = gama * pi / 180 ;
%          t_axis = Distance * sin ( gamaRadian ) + 0.5 * size ( Radon , 1 ) + 1 ; 
%          theta_axis = round ( beta + gama ) ; 
%          if ( theta_axis > 360 )
%                 theta_axis = theta_axis - 360 ; 
%          elseif ( theta_axis <= 0 )
%                 theta_axis = theta_axis + 360 ;  
%          end
%          R ( i , j ) = Radon ( round ( t_axis ) , theta_axis ) ;
%          R ( i , j ) = R ( i , j ) * Distance * cos ( gamaRadian ) ;
%      end
% end
% figure , imagesc( R ) ; 

%% direct fan

%step one 

% for i = 1 :  length ( BetaScanRange ) 
%     beta =  BetaScanRange ( i ) ;
%     betaRadian = beta * pi / 180 ;
%     source_i = Center_i + Distance * cos ( betaRadian + pi / 2 ) ;  source_j = Center_j + Distance * sin ( betaRadian + pi / 2 ) ;  %define the source in matlab coordinate
%     for j = 1 : length ( GamaDetectRange )
%                 gama =  GamaDetectRange ( j ) ; 
%                 gamaRadian = gama * pi / 180 ;
%                 Smax = 2 * Distance ; 
%                 DetectPoint_iend = Center_i + Smax * sin ( gamaRadian + betaRadian ) - Distance * sin ( betaRadian ) ;   % define end detect point in matlab coordinate
%                 DetectPoint_jend = Center_j - Smax * cos ( gamaRadian + betaRadian ) + Distance * cos ( betaRadian ) ;
%                 
%                 if ( floor ( DetectPoint_iend ) == floor ( source_i ) || ceil ( DetectPoint_iend ) == ceil ( source_i ) )                     % define projection line direction and range of i and j
%                      i_range = [ ] ;
%                 elseif ( DetectPoint_iend > source_i )                                                         
%                      i_signal = 1 ;   
%                      i_range = ceil ( source_i ) : floor ( DetectPoint_iend ) ;
%                 elseif ( DetectPoint_iend < source_i )
%                     i_signal = -1 ; 
%                     i_range = ceil ( DetectPoint_iend ) : floor ( source_i ) ;
%                 end
%                 
%                 if ( floor ( DetectPoint_jend ) == floor ( source_j ) || ceil ( DetectPoint_jend ) == ceil ( source_j ) )
%                      j_range = [ ] ;
%                 elseif ( DetectPoint_jend > source_j )
%                     j_signal = 1 ; 
%                     j_range = ceil ( source_j ) : floor ( DetectPoint_jend ) ;
%                 elseif ( DetectPoint_jend < source_j )
%                     j_signal = -1 ; 
%                     j_range = ceil ( DetectPoint_jend ) : floor ( source_j ) ;
%                 end
%                 
%                 Li = length ( i_range ) ; Lj = length ( j_range ) ;
%                 ProjectionLine = zeros ( 3 , Li  + Lj ) ;
%                 ProjectionLine ( 1 , 1 : Li  ) = abs ( ( i_range - source_i ) / ( DetectPoint_iend - source_i ) ) * Smax ;
%                 ProjectionLine ( 2 , 1 : Li  ) = i_range ;
%                 ProjectionLine ( 3 , 1 : Li  ) = 1 ;                              % 1 represents i
%                 ProjectionLine ( 1 , ( Li  + 1 ) : ( Li  + Lj ) ) = abs ( ( j_range - source_j ) / ( DetectPoint_jend - source_j ) ) * Smax ;
%                 ProjectionLine ( 2 , ( Li  + 1 ) : ( Li  + Lj ) ) = j_range ;
%                 ProjectionLine ( 3 , ( Li  + 1 ) : ( Li  + Lj ) ) = 2 ;                              % 2 represents j
%                 ProjectionLine = sortrows ( ProjectionLine' )' ;
%                 DetectPoint_i = ceil ( source_i ) ; DetectPoint_j = ceil ( source_j ) ;                   % start point of the line        
%                 
%                 for n = 1 : length ( ProjectionLine )
%                     
%                     if ( ProjectionLine ( 3 , n ) == 1 && Li ~= 0 )                                                              % define the pixel 
%                             DetectPoint_i = ProjectionLine ( 2 , n ) + 0.5 + i_signal * 0.5 ; 
%                     elseif ( ProjectionLine ( 3 , n ) == 2 && Lj ~= 0 )
%                             DetectPoint_j = ProjectionLine ( 2 , n ) + 0.5 + j_signal * 0.5 ; 
%                     end
%                     
%                     if ( Li == 0 && source_i == 0 && DetectPoint_j > 0 && DetectPoint_j <= width )                % boundary condition
%                         f = pic ( 1 , DetectPoint_j ) / 2 ;
%                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                         R ( i , j ) = R ( i , j ) + f * LengthUnit ;    
%                     end   
%                     if ( Lj == 0 && source_j == 0 && DetectPoint_i > 0 && DetectPoint_i <= width )
%                         f = pic ( DetectPoint_i , 1 ) / 2 ;
%                          LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                         R ( i , j ) = R ( i , j ) + f * LengthUnit ;    
%                     end   
%                     
%                     if ( DetectPoint_i > 0 && DetectPoint_i <= height && DetectPoint_j > 0 && DetectPoint_j <= width) 
%                          
%                          if ( Li == 0 && source_i == DetectPoint_i )
%                              if ( source_i == height )
%                                  f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
%                              else 
%                                  f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i + 1 , DetectPoint_j ) ) / 2 ;

%                              end
%                          elseif ( Lj == 0 && source_j == DetectPoint_j ) 
%                              if ( source_j == width )
%                                  f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
%                              else 
%                                  f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i , DetectPoint_j + 1 ) ) / 2 ;
%                              end
%                          else 
%                                 f = pic ( DetectPoint_i , DetectPoint_j ) ;
%                          end
% 
%                          LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
%                      end 
%                      R ( i , j ) = R ( i , j ) + f * LengthUnit ;    
%                end                
%                R ( i , j ) = R ( i , j ) * Distance * cos ( gamaRadian ) ;
%     end
% end
% R = R .* ( Distance * repmat ( cos ( GamaDetectRange )' , LBeta , 1 ) ) ;    % Pre-Weighting
% figure , imagesc( R ) ; 
            for i = 1 : LBeta                    
                    beta = BetaScanRange ( i ) ;
                    for j = 1 : LGama
                            gama = GamaDetectRange ( j ) ;

                            W = ParkFunction ( beta , StartAngle , gama , MaxGama ) ;
                            
                            R ( i , j ) = R ( i , j ) * Distance * W * cos ( gama ) ;
                            
                    end
            end

%% filter convolution based R-L
%%step two

% for i = 1 : LGama
%     N =  LGama - i ;
%         if N == 0
%             G ( i ) = 1 / ( 4 * GamaDetectInt ^2 ) ;    % radian 
%         elseif rem ( N , 2 ) == 0          %even
%             G ( LGama + N ) = 0 ;
%             G ( LGama - N ) = 0 ;
%         else                                         %odd
%            G ( LGama + N ) = - 1 / ( pi * sin ( N * GamaDetectInt )  )  ^ 2 ;
%            G ( LGama - N ) = - 1 / ( pi * sin ( N * GamaDetectInt )  )  ^ 2 ;
%         end
% end

% Rcov = zeros ( LBeta ,  LGama ) ;
% for i = 1 :  LBeta
%         cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = GamaDetectInt * cov ( LGama - 4 : 2 * LGama - 5 ) ;
% end

%Hamming window to filter

% Hamming = zeros ( length ( GamaDetectRange )   , 1 ) ; 
% Hamming ( length ( GamaDetectRange ) / 2 - 4 : length ( GamaDetectRange ) / 2 + 5 ) = hamming ( 10 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( Rcov ( i , : ) , Hamming' ) ;                              % convolution with filter
%         Rcov ( i , : ) = cov ( round ( length ( GamaDetectRange ) * 0.5 ) : round ( length ( GamaDetectRange ) * 0.5 ) + length ( GamaDetectRange ) - 1 ) / Hammingsum ;
% end

% figure , imagesc( Rcov ) ; 

%% filter convolution based S-L
%%step two

for i = 1 : LGama
    N =  LGama - i ;
        if N == 0
            G ( i ) = 2 / ( pi * GamaDetectInt )^2  ;    % radian 
        else                                         %odd
           G ( LGama + N ) = - 2 / ( pi * GamaDetectInt )^ 2 / ( 4 * N^2 - 1 ) ;
           G ( LGama - N ) = - 2 / ( pi * GamaDetectInt )^ 2 / ( 4 * N^2 - 1 ) ;
        end
end

Rcov = zeros ( LBeta ,  LGama ) ;
for i = 1 :  LBeta
        cov = conv ( R ( i , : ) , G' ) ;                              % convolution with filter
%         Rcov ( i , : ) = GamaDetectInt * cov ( LGama - 5 : 2 * LGama - 6 ) ;
        Rcov ( i , : ) = GamaDetectInt * cov ( LGama  : 2 * LGama - 1 ) ;
end
% figure , imagesc( Rcov ) ; 

%% reconstruction
%step three
% 
for i = 1 :  height
    for j = 1 : width
            RecPoint_i = i - 0.5 - Center_i ; RecPoint_j = j - 0.5 - Center_j ;                                % reconstructed point in matlab cordinate
            for betas = 1 : LBeta
                    beta =  BetaScanRange ( betas ) ;
                    t =  RecPoint_j * cos ( beta ) + RecPoint_i * sin ( beta ) ;
                    s =  -RecPoint_j * sin ( beta ) + RecPoint_i * cos ( beta ) ; 
                    s_shift = s - Distance ;
                    gama_b = atan2 ( t , abs ( s_shift ) ) ;
                    L2 = t^2 + s_shift^2 ;
                    
                    if ( gama_b >= min ( GamaDetectRange ) && gama_b < max ( GamaDetectRange ) )
                        
                            gama1index = floor (  abs ( gama_b - GamaDetectRange ( 1 ) )  / GamaDetectInt )  + 1 ;
                            gama2index = gama1index + 1 ;
                            gama1 = GamaDetectRange ( gama1index ) ; 
                            gama2 = GamaDetectRange ( gama2index ) ;
                            
                            Projection ( i , j ) = Projection ( i , j ) +  ( abs ( gama2 - gama_b ) * Rcov ( betas , gama1index ) + abs ( gama_b - gama1 ) * Rcov ( betas , gama2index ) ) / ( L2 * GamaDetectInt );       % bilinear interpolation                 
%                             Projection ( i , j ) = Projection ( i , j ) +  Rcov ( betas , ( LGama - gama1index  + 1 ) ) / ( L2 ) ;  
%                             Projection ( i , j ) = Projection ( i , j ) +  Rcov ( betas , gama1index ) / ( L2 ) ;  
                    end
            end
            Projection ( i , j ) = Projection ( i , j ) * BetaScanInt ;
    end
end
figure , imshow( fliplr ( Projection )' , [0,2000] ) ; 
title( ' Reconstructed Image ' ) ;

% figure , plot ( 1 : height , Projection ( height / 2 , : ) , 1 : height , pic ( height / 2 , : ) ) ;
% title('grey distrubition') ;

%% evaluation function
% aver=sum ( sum ( pic ) ) / ( size ( pic , 1 ) * size ( pic , 2 ) ) ;
% pd=double ( pic ) ;
% d= ( sum ( sum ( ( pd - Projection ).^2 ) ) / sum ( sum ( ( pd - aver ).^2 ) ) ) ^0.5 ;
toc
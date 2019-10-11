% modified on 2017/03/16 ZXZ
% modified on 2017/04/7 
tic
%% clear all;
% close all;
%%
% parameter define

pic = StandardPhantom ( "shepp_logan" , 257 ) ;
% load 'E:\ZXZ\Data\trial2D'
% pic = trial2D ; 

Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_i = Size ( 1 ) / 2 ;  Center_j = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.5 ;
t_range =  -tmax : t_int : tmax ;
Lt = length ( t_range ) ;

thetaint = deg2rad(1) ;     % theta unit 
Maxtheta = deg2rad(360) ;      
thetaRange = thetaint : thetaint : Maxtheta ;                                % radon scanning range, radian
Ltheta = length ( thetaRange ) ; 

%% my radon

R = zeros ( Lt ,  Ltheta ) ;   % create space to store fan projection
Smax = 2 * tmax ;
for theta_axis = 1 : Ltheta
%         theta_axis = 90 ;
        
        thetaradian = thetaRange ( theta_axis ) ;
        Htransform ( 1 , 1 ) = cos ( thetaradian ) ;
        Htransform ( 1 , 2 ) = -sin ( thetaradian ) ;
        Htransform ( 2 , 1 ) = sin ( thetaradian ) ;
        Htransform ( 2 , 2 ) = cos ( thetaradian ) ;
        for t_axis = 1 : Lt
%                 t_axis = 109 ; 
                T = t_range ( t_axis ) ;
                Pstart = [ -tmax ; T ] ;
                Pend = [ tmax ; T ] ; 
                Resultstart = Htransform * Pstart ;
                Resultend = Htransform * Pend ;
                dect_istart = Resultstart ( 1 ) + Center_i ;
                dect_jstart = Resultstart ( 2 ) + Center_j ;
                dect_iend = Resultend ( 1 ) + Center_i ;
                dect_jend = Resultend ( 2 ) + Center_j ;
                
                % define projection line direction and range of i and j
                if ( floor ( dect_iend ) == floor ( dect_istart ) || ceil ( dect_iend ) == ceil ( dect_istart ) )    
                    i_range = [ ] ;
                elseif ( dect_iend > dect_istart )                                                         
                     i_signal = 1 ;   
                     i_range = ceil ( dect_istart ) : floor ( dect_iend ) ;
                elseif ( dect_iend < dect_istart )
                    i_signal = -1 ; 
                    i_range = ceil ( dect_iend ) : floor ( dect_istart ) ;
                end
                if ( floor ( dect_jend ) == floor ( dect_jstart ) || ceil ( dect_jend ) == ceil ( dect_jstart ) )
                     j_range = [ ] ;
                elseif ( dect_jend > dect_jstart )
                    j_signal = 1 ; 
                    j_range = ceil ( dect_jstart ) : floor ( dect_jend ) ;
                elseif ( dect_jend < dect_jstart )
                    j_signal = -1 ; 
                    j_range = ceil ( dect_jend ) : floor ( dect_jstart ) ;
                end       
                
                Li = length ( i_range ) ; Lj = length ( j_range ) ;
                ProjectionLine = zeros ( 3 , Li + Lj ) ;
                ProjectionLine ( 1 , 1 : Li ) = abs ( ( i_range - dect_istart ) / ( dect_iend - dect_istart ) ) * Smax ;
                ProjectionLine ( 2 , 1 : Li ) = i_range ;
                ProjectionLine ( 3 , 1 : Li ) = 1 ;                              % 1 represents i
                ProjectionLine ( 1 , ( Li + 1 ) : ( Li + Lj ) ) = abs ( ( j_range - dect_jstart ) / ( dect_jend - dect_jstart ) ) * Smax ;
                ProjectionLine ( 2 , ( Li + 1 ) : ( Li + Lj ) ) = j_range ;
                ProjectionLine ( 3 , ( Li + 1 ) : ( Li + Lj ) ) = 2 ;                              % 2 represents j
                ProjectionLine = sortrows ( ProjectionLine' )' ;
                DetectPoint_i = ceil ( dect_istart ) ; DetectPoint_j = ceil ( dect_jstart ) ;                   % start point of the line        

                for n = 1 : length ( ProjectionLine ) - 1                                                                                    
                    if ( ProjectionLine ( 3 , n ) == 1 && Li ~= 0 )                                                               % define the pixel 
                            DetectPoint_i = ProjectionLine ( 2 , n ) + 0.5 + i_signal * 0.5 ; 
                    elseif ( ProjectionLine ( 3 , n ) == 2 && Lj ~= 0 )
                            DetectPoint_j = ProjectionLine ( 2 , n ) + 0.5 + j_signal * 0.5 ; 
                    end
                    
                    if ( Li == 0 && dect_istart == 0 && DetectPoint_j > 0 && DetectPoint_j <= width )                % boundary condition
                        f = pic ( 1 , DetectPoint_j ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end   
                    if ( Lj == 0 && dect_jstart == 0 && DetectPoint_i > 0 && DetectPoint_i <= width )
                        f = pic ( DetectPoint_i , 1 ) / 2 ;
                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end   
                    
                    if ( DetectPoint_i > 0 && DetectPoint_i <= height && DetectPoint_j > 0 && DetectPoint_j <= width ) 
                         
                        if ( Li == 0 && dect_istart == DetectPoint_i )                        % condition when projection line overlay the axis line                
                             if ( dect_istart == height )
                                 f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i + 1 , DetectPoint_j ) ) / 2 ;
                             end
                         elseif ( Lj == 0 && dect_jstart == DetectPoint_j ) 
                             if ( dect_jstart == width )
                                 f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i , DetectPoint_j + 1 ) ) / 2 ;
                             end
                         else 
                                f = pic ( DetectPoint_i , DetectPoint_j ) ;
                         end
                       
                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;                     
                         R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end
                end
        end
end
figure , imshow ( R , [] ) ;             

%%    raodn derived from fan projection 
% BetaScanInt = 1 ;             % scanning internal
% GamaDetectInt = 0.1 ;              % detective internal
% MaxGama = 45 ;                     
% MaxBeta = 360 ; 
% BetaScanRange = 0 : BetaScanInt : MaxBeta ;     % scanning range , angle between SO and aixs Y
% GamaDetectRange = -MaxGama : GamaDetectInt : MaxGama ;    % detective range
% Rfan = zeros ( length ( BetaScanRange ) ,  length ( GamaDetectRange ) ) ;   % create space to store fan projection
% Distance = max ( size ( p ) ) * 0.75 ;
% distanceint = 1 ;     % when disint = 1 , total time =  103.665 s ; 
% for i = 1 :  length ( BetaScanRange )
%     beta =  BetaScanRange ( i ) ;
%     betaRadian = beta * pi / 180 ;
%     dect_istart = Center_i - Distance * cos ( betaRadian ) ;  dect_jstart = Center_j - Distance * sin ( betaRadian ) ;  %deine the source
%     for j = 1 : length ( GamaDetectRange )
%                 gama =  GamaDetectRange ( j ) ; 
%                 gamaRadian = gama * pi / 180 ;
%                 S = 0 ;              % initialize the raylength 
%                while ( S < 2 * Distance )
%                         CO = sqrt ( ( S * cos ( gamaRadian ) - Distance) ^ 2 +  ( S * sin ( gamaRadian ) ) ^2 ) ;                   
%                             DetectPoint_i = Center_i - S * sin ( gamaRadian + betaRadian - pi / 2 ) - Distance * sin ( betaRadian + pi / 2 ) ;   % define detect point
%                             DetectPoint_j = Center_j + S * cos ( gamaRadian + betaRadian - pi / 2 ) + Distance * cos ( betaRadian + pi / 2 ) ;
%                         if ( DetectPoint_i > 1 && DetectPoint_i < height && DetectPoint_j > 1 && DetectPoint_j < width )      % ensure the projection from the object
% %                             iLU = floor ( DetectPoint_i ) ; jLU = floor ( DetectPoint_j ) ;                                            
% %                             iLD = floor ( DetectPoint_i ) +1 ; jLD = floor ( DetectPoint_j ) ;
% %                             iRU = floor ( DetectPoint_i ) ; jRU = floor ( DetectPoint_j ) + 1 ;                                            
% %                             iRD = floor ( DetectPoint_i ) +1 ; jRD = floor ( DetectPoint_j ) + 1 ;
% %                             fU = ( jRU - DetectPoint_j ) * pic ( iLU , jLU ) +  ( DetectPoint_j - jLU ) * pic ( iRU , jRU ) ;
% %                             fD = ( jRD - DetectPoint_j ) * pic ( iLD , jLD ) +  ( DetectPoint_j - jLD ) * pic ( iRD , jRD ) ;
% %                             f = (  iLD - DetectPoint_i  ) * fU + ( DetectPoint_i - iLU ) * fD ;                             % bilinear interpolation 
%                                 DetectPoint_i = round ( DetectPoint_i ) ; DetectPoint_j = round ( DetectPoint_j ) ; 
%                                 f = p ( DetectPoint_i , DetectPoint_j ) ;
%                             Rfan ( i , j ) = Rfan ( i , j ) + f * distanceint ;
%                         end 
%                         S = S + distanceint ;
%                end 
%     end
% end
% figure , imagesc( Rfan ) ; 
% 
% tmax = round ( 0.5 * sqrt ( 2 ) * size ( p , 1 ) ) ;
% R = zeros ( 2 * tmax + 1 ,  length ( theta ) ) ;   % create space to store fan projection
% 
% for theta_axis = theta
%     thetaradian = theta_axis * pi / 180 ;
%     for t_axis = 1 : 2 * tmax + 1
%         gamaradon = asin ( ( t_axis - tmax -1 ) / Distance );
%         Betaradon = round ( ( thetaradian - gamaradon ) * 180 / pi ) ;
%         if ( Betaradon > 360 )
%             Betaradon = Betaradon -360 ;
%         elseif ( Betaradon < 1)
%             Betaradon = Betaradon + 360 ;
%         end
%         gamaradon = round ( gamaradon * 180 / pi / GamaDetectInt )+ ( length ( GamaDetectRange) + 1 ) / 2 ; 
%         if ( gamaradon >= 1 && ( gamaradon <= length ( GamaDetectRange ) ) && Betaradon >= 1 && ( Betaradon <= length ( BetaScanRange ) ) )
%             R ( t_axis , theta_axis ) = Rfan ( Betaradon , gamaradon ) ; 
%         end
%     end
% end
% figure,imagesc(R);        

%%  inner radon
% [ R , xp ] = radon ( pic , thetaRange ) ;
% figure , imagesc ( R ) ;                               
% title ( ' Radon transform image ' ) ;
%% formualr
%         A    a     b    x0    y0    phi
%        ---------------------------------
% shep = [  2   .69   .92    0     0     0                                                      % in ground coordinate
%         -.98 .6624 .8740   0  -.0184   0
%         -.02 .3100 .1100   .22    0    72
%          -.02 .4100 .1600  -.22    0     108
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
% parfor num = 1 : Ltheta * Lt
%     
%         thetaindex = mod ( num - 1 ,  Ltheta ) + 1 ;
%         tindex = floor ( ( num - 1 ) / Ltheta ) + 1 ;
%         
%         thetaRadian =  thetaRange ( thetaindex ) * pi / 180 ;       % angle to radian            
%         t = t_range ( tindex ) ;      
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
% end
% R = reshape ( R , Ltheta , Lt ) ;    

% figure , imshow ( R , [] ) ; 

%% filter

h = zeros ( ( Lt * 2 - 1 ) , 1 ) ;
% RL
% for i = 0 : Lt - 1                %假设采样间隔为1，有限带宽斜变滤波器 
%     if i == 0
%         h ( Lt ) = 1 / ( 2 * t_int )^2 ;
%     elseif rem ( i , 2 ) == 0          % even
%         h ( Lt - i ) = 0 ;
%         h ( Lt + i ) = 0 ;
%     else
%         h ( Lt - i )= -1 / ( i * pi * t_int ) ^ 2 ;     % odd
%         h ( Lt + i ) = -1 / ( i * pi * t_int ) ^ 2 ;
%     end
% end
% SL 
for i = 0 : Lt - 1                %假设采样间隔为1，有限带宽斜变滤波器 
    if i == 0
        h ( Lt ) = 1 / ( pi * t_int )^2 ;
    else
        h ( Lt - i )= -1 / ( pi * t_int ) ^ 2 / ( 4 * i^2 - 1 ) ;     % odd
        h ( Lt + i ) = -1 / ( pi * t_int ) ^ 2 / ( 4 * i^2 - 1 ) ;
    end
end

x = zeros ( Ltheta , Lt ) ;
parfor i = 1 : Ltheta
    s = R ( : , i ) ;
    xx = t_int * conv ( s , h' ) ;         %与滤波器卷积
    x ( i , : ) = xx ( Lt : 2 * Lt - 1 ) ;     %取卷积结果中心部分
end

% Hamming = zeros ( ( H * 2 - 1 )  , 1 ) ; 
% Hamming ( H - 9 : length ( Sdomain ) / 2 + 10 ) = hamming ( 15 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% for i = 1 :  1 : N
%         cov = conv ( x ( : , i ) , Hamming' ) ;                              % convolution with filter
%         xx ( : , i ) = cov ( H : 2 * H - 1 ) / Hammingsum ;
% end

% figure,imshow(x,[]);
% title('Filtered image');
%
%% reconstruct

f = zeros ( size ( pic ) ) ;                        % store back-projection
Resolution2 =  Size(1) / height ;                                 % define the resolution of reconstruction

ProjectionShift = zeros ( Ltheta , height , width ) ;
ProjectionShift = reshape ( ProjectionShift , 1 , [] ) ;

parfor num = 1 : Ltheta * height * width                       %对每个像素点都进行反投影回抹重建
            
            i = floor ( ( num - 1 ) / ( width * Ltheta ) ) + 1 ;                % piexl number
            j = floor ( ( num - ( i - 1 ) * width * Ltheta - 1 )  / Ltheta ) + 1 ; 
            k = mod ( num - 1 , Ltheta ) + 1 ;
                       
            thetaradian = thetaRange ( k ) ;
            Point_i = ( i - 0.5 ) * Resolution2 ; Point_j = ( j - 0.5 ) * Resolution2 ;            % center point of the pixel
            
            t = -( Point_i - Center_i ) * sin ( thetaradian ) + ( Point_j - Center_j ) * cos ( thetaradian ) ;
            t1index = floor ( ( t + tmax ) / t_int ) + 1 ;
            t2index = t1index + 1 ;
            if ( t >= min ( t_range ) && t < max ( t_range ) ) 
                ProjectionShift ( num ) = ( ( t_range ( t2index ) - t ) * x ( k , t1index ) + ( t - t_range ( t1index ) ) * x ( k , t2index ) ) / t_int * thetaint ;  % interpolation
            end
end

for i = 1 : height
    for j = 1 : width
            f ( i , j )= sum ( ProjectionShift ( 1 + ( i - 1 ) * width * Ltheta + ( j - 1 ) * Ltheta : Ltheta + ( i - 1 ) * width * Ltheta + ( j - 1 ) * Ltheta ) ) ;
    end    
end       
figure , imshow ( flipud(f) , [ 0 , 0.5 ] ) ;
%figure,imshow(f,[]);
title('Reconstructed image');
%%   display

% figure,plot( 1 : size ( pic , 1 ) , f ( 257 , : ) , 1 : size ( pic , 1 ) , pic ( 257 , : ) ) ;
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 512 1 1.05 ] ) ;

% aver=sum(sum(p))/(size(p,1)*size(p,2));
% pd=double(p);
% d=(sum(sum((pd-f).^2))/sum(sum((pd-aver).^2)))^0.5;
toc ;
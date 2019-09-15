% 2017/04/24  ZXZ
% G-FDK of saddle line
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
clear all ;
Size = [ 60 , 60 , 60 ] ;     % actual range
% pic = phantom3d ( 513 ) ;     % original picture  


PicSize = 512 ;
% [ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
t_length = PicSize ; s_length = PicSize ; z_length = PicSize ;              % store the size of picture
Resolution = max ( Size ) / t_length ;                           
Rpic = max ( Size ) * sqrt ( 3 ) / 2 ;                                         % radius of project
Rplane = max ( Size ) / 2 ;                    % radius of project in the plane
HeightTraj = 5 ;              % hieght of saddle line 

Projection = zeros ( PicSize , PicSize , PicSize ) ;      % reconstruct the pic
Resolution2 = max ( Size ) / PicSize ; 

PInt = 0.4 ;                                    % interval of P ( 0.1 exact )
MaxP = Rplane * 1.1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
LP = length ( Pdomain ) ;
 
% DisRatio = 4 ; 
Distance = 103.7 ; % Rpic * DisRatio ;                     % distance between source and center point 
 
XigamaInt = 0.4 ;                                      % interval of Xigama ( 0.1 exact )
% MaxXigama = Rplane * 1.1 ;
 MaxXigama =  Rplane * 1.1 ;
% MaxXigama = ( max ( Size ) + 2 * HeightTraj ) * Distance / ( 2 * Distance - sqrt ( 2 ) * max ( Size ) ) ; 
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
LXigama = length ( Xigamadomain ) ;
 
Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center

% FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = 1 ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = 360 ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ; 

R = zeros ( LBeta , LXigama , LP ) ;   % create space to store fan projection

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
% Disk phantom
 %         A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
% Shepp =    [  1    .7  .7   .06      0       0       0        0      
%                     1    .7  .7   .06      0       0      .24        0
%                     1    .7  .7   .06      0       0       -.24        0
%                     1    .7  .7   .06      0       0       .48        0
%                     1    .7  .7   .06      0       0       -.48        0
%                     1    .7  .7   .06      0       0       .72        0
%                     1    .7  .7   .06      0       0       -.72      0 ] ;
        
        
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
%             source_t = Center_t + Distance * cos ( betaRadian + pi / 2 ) ;      % define the source in matlab coordinate
%             source_s = Center_s + Distance * sin ( betaRadian + pi / 2 ) ;    % source locates on the S axis
%             zextra = -HeightTraj * cos ( 2 * ( betaRadian + pi / 4 ) ) ; 
%             source_z  = Center_z + zextra ;    
%             
%                     gama = atan ( Xigamadomain ( Xigamaindex ) / Distance ) ;          % radian angle in s-z coordinate plane ( is converted compared with book )
%                     r = Xigamadomain ( Xigamaindex ) * Distance / sqrt ( Distance^2 + Xigamadomain ( Xigamaindex )^2 ) ;                % parameter in parallel ray
%                     Distance_shift = Distance / cos ( -gama ) ;                    % length of DO'
%                              
%                              FanA = atan ( Pdomain ( Pindex ) / Distance ) ; 
%                              thetaRadian = betaRadian + FanA ;                                    % radian angle in s'-t coordinate plane 
%                              t = Pdomain ( Pindex ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 ) ;                % parameter in parallel ray
%                              for n = 1 : 12
%                                         A = Shepp ( n , 2 ) * Proportion ; B = Shepp ( n , 3 ) * Proportion ; C = Shepp ( n , 4 ) * Proportion ;                    % information of oval
%                                         rou = Shepp ( n , 1 ) ;
%                                         X1 = Shepp ( n , 5 ) * Proportion ; Y1 = Shepp ( n , 6 ) * Proportion ; Z1 = Shepp ( n , 7 ) * Proportion - zextra ;
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
%                              R ( num ) = R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% %                              if ( betaRadian <= 2 * ( FanAmax - FanA ) )
% %                                     W = sin ( pi * betaRadian / ( FanAmax - FanA ) / 4 )^2 ;
% %                                     R ( num ) = W * R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% %                              elseif ( betaRadian > 2 * ( FanAmax - FanA ) && betaRadian <= pi - 2 * FanA )
% %                                     R ( num ) = R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% %                              elseif ( betaRadian > pi - 2 * FanA && betaRadian <= pi + 2 * FanAmax )        
% %                                     W = sin ( pi * ( pi + 2 * FanAmax - betaRadian ) / ( FanAmax + FanA ) / 4 )^2 
% %                                     R ( num ) = W * R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( Pindex )^2 + Xigamadomain ( Xigamaindex )^2 ) ;
% %                              end  
% end
% R = reshape ( R , LBeta , LXigama , LP ) ;

% figure , imshow ( squeeze ( R ( : , 331  , : ) ) , [] ) ; 

 %% filter    S-L
 load ( ' D:\TestCpp\CT\Data\FDK\RFBPimage3\SaddleR103.9.mat ' ) ;
 G = zeros ( 2 * LP - 1  , 1) ;                                                      % create filter

for Pindex = 0 : LP - 1 
        if Pindex == 0
            G ( LP ) = 1 / ( pi * PInt )^2 ;    % radian 
        else                                       
           G ( LP + Pindex ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * Pindex^2 - 1 ) ;
           G ( LP - Pindex ) = - 1 / ( pi * PInt ) ^ 2 / ( 4 * Pindex^2 - 1 ) ;
        end
end

Rcov = zeros ( LBeta ,  LXigama , LP ) ;
parfor Betaindex = 1 :  LBeta
    for Xigamaindex = 1 : LXigama
    
        cov = conv ( squeeze ( R ( Betaindex , Xigamaindex , : ) ) , G' ) ;                              % convolution with filter
        Rcov ( Betaindex , Xigamaindex , : ) = PInt * cov ( LP : 2 * LP - 1 ) ;

    end
end

clear R G ;

% Hamming = zeros ( length ( Pdomain )  , 1 ) ; 
% Hamming ( length ( Pdomain ) / 2 - 7 : length ( Pdomain ) / 2 + 7 ) = hamming ( 15 ) ;                   % convolve with hamming window
% Hammingsum = sum ( Hamming ) ;             % to divide the sum of hamming window 
% for i = 1 :  length ( BetaScanRange )
%         cov = conv ( squeeze ( R ( i , j , : ) ) , Hamming' ) ;                              % convolution with filter
%         Rcov ( i , j , : ) = cov ( round ( length ( Pdomain ) * 0.5 ) : round ( length ( Pdomain ) * 0.5 ) + length ( Pdomain ) - 1 ) / Hammingsum ;
% end

% figure , imshow ( squeeze ( Rcov ( : , 331 , : ) ) , [] ) ; 

%% reconstruct

Display = zeros ( PicSize , PicSize ) ;
ProjectionSlice = zeros ( PicSize , PicSize , LBeta ) ;
ProjectionSlice = reshape ( ProjectionSlice , 1 , [] ) ;
% for t = 1 : t_length
    
    t = 257 ;
    
    parfor num = 1 : t_length * z_length * LBeta
%                 s = 257 ;
               z = floor ( ( num - 1 ) / ( LBeta * t_length ) ) + 1 ;
%                z = 509 ;
               s = floor ( ( num - ( z - 1 ) * t_length * LBeta - 1 ) / LBeta ) + 1 ; 
%                     betas = 1 ;
               betas = mod ( num - 1 , LBeta ) + 1 ;
               
                    beta =  BetaScanRange ( betas ) ;
                    betaRadian = beta * pi / 180 ;
                    
                    source_t = Center_t + Distance * cos ( betaRadian + pi / 2 ) ;                           % define the source
                    source_s = Center_s + Distance * sin ( betaRadian + pi / 2 ) ;   
                    zextra = -HeightTraj * cos ( 2 * ( betaRadian + pi / 4 ) ) ; 
                    source_z  = Center_z + zextra ;    
                    
                    image_t = ( t - 0.5 ) * Resolution2 - Center_t  ;  image_s = ( s - 0.5 ) * Resolution2 - Center_s  ; image_z = ( z - 0.5 ) * Resolution2 - Center_z  ;           % image pixel in ground coordinate
                    
                    dect_t = image_t * cos ( betaRadian ) + image_s * sin ( betaRadian ) ;          % in rotate coordinate
                    dect_s = - image_t * sin ( betaRadian ) + image_s * cos ( betaRadian ) ; 
                    dect_z = image_z ;     
                    
                    LengthRatio = Distance / ( Distance - dect_s ) ; 
                    Xigama_domain = ( dect_z - zextra ) * LengthRatio ; 
                    P_domain = dect_t * LengthRatio ; 
                    
                    Xigama_domainN1 =  floor ( ( Xigama_domain + MaxXigama ) / XigamaInt ) + 1 ; 
                    Xigama_domainN2 = Xigama_domainN1 + 1 ;
                    P_domainN1 =  floor ( ( P_domain + MaxP ) / PInt ) + 1 ; 
                    P_domainN2 = P_domainN1 + 1 ;

                    if  ( P_domain >= -MaxP && P_domain < max ( Pdomain )  && Xigama_domain >= -MaxXigama && Xigama_domain < max ( Xigamadomain ) )
                            P_domain1 = Pdomain ( P_domainN1 ) ; P_domain2 = Pdomain ( P_domainN2 ) ;
                            Xigama_domain1 = Xigamadomain ( Xigama_domainN1 ) ; Xigama_domain2 = Xigamadomain ( Xigama_domainN2 ) ;  
                            
                            % bilinear interpolation
                            Xig1 = Xigama_domain - Xigama_domain1 ; Xig2 = Xigama_domain2 - Xigama_domain ;
                            P1 = P_domain - P_domain1 ; P2 = P_domain2 - P_domain ;
                            ProjectionSlice ( num ) = ( Xig2 * P2 * Rcov ( betas , Xigama_domainN1 , P_domainN1 ) ...
                                + Xig1 * P2 * Rcov ( betas , Xigama_domainN2 , P_domainN1 ) + Xig2 * P1 * Rcov ( betas , Xigama_domainN1 , P_domainN2 ) ...
                                + Xig1 * P1 * Rcov ( betas , Xigama_domainN2 , P_domainN2 ) ) / ( PInt * XigamaInt ) * LengthRatio^2 * BetaScanInt * pi / 180 ;
                    end                       
    end
     
% end
for  z = 1 : z_length 
    for s = 1 : t_length 
        Display ( z , s ) = sum ( ProjectionSlice ( 1 + ( ( z - 1 ) * t_length * LBeta + ( s - 1 ) * LBeta ) : LBeta + ( ( z - 1 ) * t_length * LBeta + ( s - 1 ) * LBeta ) ) ) ;
    end
end   
clear Rcov ProjectionSlice ;
% figure , imshow ( flipud ( Display ) , [  ] ) ; 
figure , imshow ( Display , [ 0.15 0.35 ] ) ; 
Time = toc ;
% figure , plot ( 1 : 513 , squeeze ( Projection (  , : , 50 ) ) , 1 : 513 , squeeze ( pic ( size ( pic , 1 ) / 2 , : , 50 ) )  ) ;
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 513 1 1.05 ] ) ;
% set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;
% set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;
%% read fig
% a = open ( 'D:\TestCpp\CT\Data\FDK\Modified_Para\illed\t257.fig' ) ;
% h = get (gca , 'Children') ;
% data = get ( h , 'Cdata') ;
% figure,imshow ( squeeze ( pic ( 257 , : , : ) )' , [1,1.05] )     s=257

% x=[1:513] ;
% y= -8.116e-7 .* ( x -257 ).^2 ;
% figure , plot ( 1 : 513  , err ( : , 125 ) , 1 : 513 , y );axis ( [ 0 513 -0.05 0 ] ) ;


% a = squeeze ( pic ( 257 , : , : ) )' ;
% figure , imshow ( flipud ( a ) , [1,1.05] ) ;















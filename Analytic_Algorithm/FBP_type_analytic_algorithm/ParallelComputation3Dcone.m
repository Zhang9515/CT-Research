% parpool  local ;
%2016/12/03 ZXZ
%2017/4/8 modified
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
clear all ;
Size = [ 60 , 60 , 60 ] ;     % actual range
% pic = phantom3d ( 513 ) ;     % original picture  

BetaScanInt = 1 ;             % scanning internal    ( 0.3 exact )          
MaxBeta = 360 ; 
BetaScanRange = BetaScanInt + 0.1 : BetaScanInt : MaxBeta + 0.1  ;     % scanning range , angle between SO and aixs Y
LBeta = length ( BetaScanRange ) ;

% [ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
t_length = 513 ; s_length = 513 ; z_length = 513 ;              % store the size of picture
Resolution = max ( Size ) / t_length ;                           
Rpic = max ( Size ) * sqrt ( 3 ) / 2 ;                                         % radius of project

Projection = zeros ( 513 , 513 , 513 ) ;      %reconstruct the pic
Resolution2 = max ( Size ) / 513 ; 

PInt = 0.5 ;                                    %interval of P ( 0.1 exact )
MaxP = Rpic ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
LP = length ( Pdomain ) ;

XigamaInt = 0.5 ;                                      % interval of Xigama ( 0.1 exact )
MaxXigama = Rpic ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

DisRatio = 4 ; 
Distance = Rpic * DisRatio ;                     % distance between source and center point

R = zeros ( LP , LXigama , LBeta ) ;   % create space to store fan projection
R = reshape ( R , 1 , [] ) ;

%%  formular projection  P.104
%                    A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
Shepp =    [        2  .6900  .920  .900      0       0       0      0    
                 -.98  .6624  .874  .880      0       0       0      0    
                 -.02  .310  .1100  .220    .22      0      -.25     72    
                -.02  .410  .1600  .280   -.22       0     -.25    108     
                 .02  .0460  .046  .046      0      .1     -.25     0     
                 .02  .2100  .250  .500      0      .35    -.25     0     
                 .01  .0460  .046  .046      0     -.1     -.25     0      
                 .01  .0460  .023  .020   -.08    -.605    -.25     90    
                 .01  .0230  .023  .020      0   -.606     -.25     90   
                 .01  .0460  .023  .020    .06   -.605     -.25      0  
                 .02  .0560  .040  .10    .06    -.105     .625      0
                -.02  .0560  .056  .10       0     .105    .625      0 ];

Proportion = max ( Size ) / 2 ;
parfor num = 1 :  LBeta * LXigama * LP 
            i = floor ( ( num -1 ) / ( LXigama * LP ) ) + 1 ; 
            beta =  BetaScanRange ( i ) ;
            betaRadian = beta * pi / 180 ;                                          % angle to radian
            source_t = Center_t - Distance * sin ( betaRadian ) ;      % define the source in ground coordinate
            source_s = Center_s + Distance * cos ( betaRadian ) ;    % source locates on the S axis
            source_z  = Center_z ;                                                        
                    j = floor ( ( num - ( i - 1 ) * LXigama * LP - 1 ) / LP ) + 1 ; 
                    gama = atan ( Xigamadomain ( j ) / Distance ) ;          % radian angle in s-z coordinate plane
                    r = Xigamadomain ( j ) * Distance / sqrt ( Distance^2 + Xigamadomain ( j )^2 ) ;                % parameter in parallel ray
                    Distance_shift = Distance / cos ( gama ) ;                    % length of DO'
                             k = mod ( num - 1 , LP ) + 1 ; 
                             thetaRadian = betaRadian + atan ( Pdomain ( k ) / Distance ) ;                                    % radian angle in s'-t coordinate plane 
                             t = Pdomain ( k ) * Distance / sqrt ( Distance^2 + Pdomain ( k )^2 ) ;                % parameter in parallel ray
                             for n = 1 : 12
                                        A = Shepp ( n , 2 ) * Proportion ; B = Shepp ( n , 3 ) * Proportion ; C = Shepp ( n , 4 ) * Proportion ;                    % information of oval
                                        rou = Shepp ( n , 1 ) ;
                                        X1 = Shepp ( n , 5 ) * Proportion ; Y1 = Shepp ( n , 6 ) * Proportion ; Z1 = Shepp ( n , 7 ) * Proportion ;
                                        fai = Shepp ( n , 8 ) * pi / 180 ;
                                        thetas = thetaRadian - fai ;
                                        ts = t - X1 * cos ( thetaRadian ) - Y1 * sin ( thetaRadian ) ;
                                        rs = r + ( Y1 * cos ( thetaRadian ) - X1 * sin ( thetaRadian ) ) * sin ( gama ) - Z1 * cos ( gama ) ;   % maybe gama minus
                                        a2 = C^2 * ( ( B * sin ( thetas ) )^2 + ( A * cos ( thetas ) )^2 ) * ( cos ( gama ) )^2 + ( A * B * sin ( gama ) )^2 ;
                                        Sq = a2 - ts^2 * ( ( C * cos ( gama ) )^2 + ( ( B * cos ( thetas ) )^2 + ( A * sin ( thetas ) )^2 ) * ( sin ( gama ) )^2 ) ...
                                                - ( rs )^2 * ( ( B * sin ( thetas ) )^2 + ( A * cos ( thetas ) )^2 ) * ( 7 + cos ( 4 * gama ) ) / 8 ...
                                                - 2 * ts * rs * sin ( gama ) * cos ( thetas ) * sin ( thetas ) * ( B^2 - A^2 ) ;
                                        if ( Sq >= 0 )
                                                R ( num ) = R ( num ) + 2 * rou * A * B * C * sqrt ( Sq ) / a2 ;
                                        end                            
                             end
                             R ( num ) = R ( num ) * Distance / sqrt ( Distance^2 + Pdomain ( k )^2 + Xigamadomain ( j )^2 ) ;
end
R = reshape ( R , LP , LXigama , LBeta ) ;
% figure , imshow ( squeeze ( R ( : , 105  , : ) ) , [] ) ; 
Time = toc ;
% delete ( gcp ( 'nocreate' ) ) ;
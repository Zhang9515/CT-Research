%2018/07/25 by ZXZ
% iterative FDK algorithim
tic 
clear;
Size = [ 60 ; 60 ; 60 ] ;     % actual range 60
pic = zeros ( 128 ) ;

pic = single(pic);

% PicSize = 512 ;
[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
% t_length = PicSize ; s_length = PicSize ; z_length = PicSize ;              % store the size of picture
Resolution = max ( Size ) / t_length ;                           
Rpic = max ( Size ) * sqrt ( 3 ) / 2 ;                                         % radius of project (51.9615 for size 60)
Rplane = max ( Size ) * sqrt ( 2 ) / 2 ;                    % radius of project in the plane

% Projection = zeros ( PicSize , PicSize , PicSize ) ;      %reconstruct the pic
Resolution2 = max ( Size ) / size ( pic,1 ) ; 

% PInt = 0.4 ;                                    %interval of P ( 0.1 exact )
PInt = 0.4;
MaxP = Rplane * 1.1 ;
% MaxP = PInt * 83 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

DisRatio = 4 ; 
Distance = 200; % Rpic * DisRatio % 207.6;                     % distance between source and center point
% Distance = 730 ; 

% FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = deg2rad(1) ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = deg2rad(360) ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 
% BetaScanRange has been radian

dU = Resolution * 1.0 ;
LU = floor(2 * Rpic / dU );
Udomain = Rpic - dU * (1 : LU);
Udomain = single(Udomain');

ErrorSlice = 0 ;

t = 64 ;

betas = 1 ;     % in ground image coordinate

z = 1 ; 

s = 1;

betaRadian =  BetaScanRange ( betas ) ;

image_t = ( t - 0.5 ) * Resolution2 - Center_t  ;  image_s = ( s - 0.5 ) * Resolution2 - Center_s  ; image_z = ( z - 0.5 ) * Resolution2 - Center_z  ;           % image pixel in ground coordinate

dect_t = image_t * cos ( betaRadian ) + image_s * sin ( betaRadian ) ;          % in rotate coordinate
dect_s = - image_t * sin ( betaRadian ) + image_s * cos ( betaRadian ) ; 
dect_z = image_z ;     

LengthRatio = Distance / ( Distance - dect_s ) ; 

        for u = 1 : LU
            for p = 1 : LP      % it finds that if first p then u, the result is wrong. if first u then p, the result is right

            T_deriv = Pdomain(p) * ( Distance - Udomain(u)) / Distance;                         % in rotate coordinate
            S_deriv = Udomain(u);

            X_deriv = T_deriv * cos ( betaRadian ) - S_deriv * sin ( betaRadian ) ;          % in ground coordinate
            Y_deriv = T_deriv * sin ( betaRadian ) + S_deriv * cos ( betaRadian ) ;

            X_deriv_index = floor( ( X_deriv + Center_t) / Resolution2)+1;   

            if  ( X_deriv_index >= 1 && X_deriv_index <= t_length  )

                ErrorSlice = ErrorSlice + round(Pdomain(p));    
%                     ErrorSlice = ErrorSlice + p ;  
                    
            end  % if- index 

        end % LP    
end % LU

disp(ErrorSlice)


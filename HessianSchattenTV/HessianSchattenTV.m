% 2018/11/14
clear all
tic
Size = [ 64 , 64 , 64 ] ;     % actual range 60

pic = phantom3d ( 'Modified Shepp-Logan' , 64 ) ;     % original picture  
% load ( ' D:\TestCpp\CT\Data\FDK\SFBPimage3.mat ' ) ;       % head ct data
% pic = double ( SFBPimage ) ;
% load ( ' D:\TestCpp\CT\Data\FDK\Chest256_160.mat ' ) ;       % chest ct data
% pic = double ( Chest ) ;

pic = single(pic);

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture
%t_length = 64 ; s_length = 64 ; z_length = 64 ;              % store the size of picture
PicSize = [ t_length , s_length , z_length ] ;

Resolution = Size ./ PicSize ;                           
Rplane = max ( Size ) / 2 ;                    % radius of project in the plane

% PInt = 0.4 ;                                    %interval of P ( 0.1 exact )
PInt = 1 ;
MaxP = Rplane * 1.1 ;
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

% XigamaInt = 0.4 ;                                      % interval of Xigama ( 0.1 exact )
XigamaInt = 1 ; 
MaxXigama = Rplane * 1.1 ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = single(Xigamadomain');
LXigama = length ( Xigamadomain ) ;

Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

DisRatio = 4 ; 
Distance = 103.9; % Rpic * DisRatio ;                     % distance between source and center point

FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = 1.5 ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = 360 ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

% CovarianceMatrix = eye( LBeta * LXigama * LP ) ;

%% 
alpha = 3e4 ; Tao = 1 ; LipschitzConstant = 144 * Tao^2 ; ps =3; Knest0 = 1; Kmm0 = 1; ProjectionOnBallofBsq0 = zeros( 3 , 3 , t_length * s_length * z_length ) ;
Tmm = 10;   % outerloop
Tnest = 10 ; % innerloop

picvector = reshape (pic, t_length * s_length * z_length, 1) ;
V = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;

Kmmprevious = Kmm0 ;
for tmm = 1 : Tmm
    ProjectionData = ProjectionCone_3D (Ucurrent, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;

    Residual = Backprojection( V - ProjectionData , Xigamadomain , Pdomain , BetaScanRange , Distance, Size, t_length, s_length, z_length) ;

    Z = Ucurrent + alpha^(-1) * Residual ; 

    ProjectionOnRN = Z - Tao * HstarOmega( ProjectionOnBallofBsq , ps) ; 

    GradientofG = Tao * HessianLOG3D( ProjectionOnRN ) ;

      %% Nesterov method: orthogonal projections , to update omega, we use the result of ProjectionOnBallofBsq as 
      % Omega
    Knestprevious = Knest0 ; ProjectionOnBallofBsqprevious = ProjectionOnBallofBsq0 ;
    for tnest = 1 : Tnest
          ProjectionOnBallofBsq = Omegaprevious + GradientofG / LipschitzConstant ; 
          for n = 1 : t_length * s_length * z_length

                % q = 2
                Fnorm = norm ( ProjectionOnBallofBsq( : , : , n ), 'fro' ) ;
                if ( Fnorm > 1 )
                    ProjectionOnBallofBsq( : , : , n ) = ProjectionOnBallofBsq( : , : , n ) / Fnorm ;
                end

                % q = inf
            %     [ Usvd, ProjectionOnBallofBsq, Vsvd ] = svd ( ProjectionOnBallofBsq( : , : , n ) ) ; 

          end %n
          Knestcurrent = ( 1 + sqrt( 1 + Knestprevious^2 ) ) / 2 ;
          Omegacurrent = ProjectionOnBallofBsqprevious + ( Knestprevious - 1 ) / Knestcurrent .* ( ProjectionOnBallofBsq - ProjectionOnBallofBsqprevious ) ;
          Knestprevious = Knestcurrent ;
          Omegaprevious = Omegacurrent;
          
    end  %nest
    S = Z - Tao * HstarOmega( OmegaCurrent , ps) ;
    Kmmcurrent = ( 1 + sqrt( 1 + Kmmprevious^2 ) ) / 2 ;
    
    
    
    Kmmprevious = Kmmcurrent ; 
end  %mm


% R = reshape ( R , LP , LXigama, LBeta ) ;
% figure,imshow3Dfull(R,[])



toc
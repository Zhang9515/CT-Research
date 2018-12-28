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
BetaScanInt = deg2rad(1) ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = deg2rad(360) ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

% CovarianceMatrix = eye( LBeta * LXigama * LP ) ;

%% FDK as initial input

picvector = reshape (pic, t_length * s_length * z_length, 1);
clear pic ;
% angle input should be radian ! 
V = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance);
V = reshape ( V , LP , LXigama, LBeta ) ;
% figure,imshow3Dfull(R,[])
% clear picvector;

V = single(V) ;
z_length = 30;
FDKresult = FDK ( V , Xigamadomain , Pdomain , BetaScanRange , Distance, Size, t_length, s_length, z_length) ;
% FDKresult = reshape( FDKresult, t_length, s_length, z_length ) ;
% figure,imshow3Dfull( FDKresult , [] )
%  
%% hybrid HS with TV iterative method

OrderofHS = 2 ; 
alpha = 2e10 ;   % 2e10 according to paper 
Tao = 1 ; LipschitzConstant = 144 * Tao^2 ; ps =3; 
Knesterov0 = 1; Kmm0 = 1; ProjectionOnBallofBsq0 = zeros( 3 , 3 , t_length * s_length * z_length ) ;
U0 = FDKresult ; 

Tmm = 40;   % outerloop
Tnesterov = 30 ; % innerloop

FDKresult_single = single(FDKresult) ;
AU = ProjectionCone_3D (FDKresult_single, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;

Kmmprevious = Kmm0 ;

FDKresult = reshape ( FDKresult , t_length , s_length , z_length ) ;
Cfai0 = ObjectiveFunction( V , AU , FDKresult, OrderofHS, Tao ) ;

Cfaiprevious = Cfai0 ; 
Uprevious = U0 ; 
Xprevious = U0 ;
outer_threshold = 1e-3 ;
inner_threshold = 1e-3 ;

for tmm = 1 : Tmm
    disp(['outerloop: ', num2str(tmm) ])
    
    Xprevious_single = single(Xprevious) ;
    ProjectionData = ProjectionCone_3D (Xprevious_single, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;
    
    ProjectionData_single = single(  V - ProjectionData ) ;
    Residual = Backprojection( ProjectionData_single , Xigamadomain , Pdomain , BetaScanRange , Distance, Size, t_length, s_length, z_length) ;
    
    qq = reshape(Residual , t_length, s_length, z_length ) ; %
    figure,imshow3Dfull(qq,[]) %
    
    Z = Xprevious + Residual / alpha ; 

    %% Nesterov method: orthogonal projections , to update omega, we use the result of ProjectionOnBallofBsq as 
      % Omega 
    Knesterovprevious = Knesterov0 ; ProjectionOnBallofBsqprevious = ProjectionOnBallofBsq0 ;
    Omegaprevious = ProjectionOnBallofBsq0 ; 
    
    for tnesterov = 1 : Tnesterov
          disp(['innerloop: ', num2str(tnesterov) ])
          
          ProjectionOnRN = Z - Tao * HstarOmega( Omegaprevious , ps) ; 
          ProjectionOnRN = reshape( ProjectionOnRN , t_length , s_length , z_length ) ;
    
          GradientofG = Tao * HessianLOG3D( ProjectionOnRN ) ;
          ProjectionOnBallofBsq = Omegaprevious + GradientofG / LipschitzConstant ; 
          for n = 1 : t_length * s_length * z_length

                % q = 2, OrderofHS = p =2 , 1/p + 1/q = 1
                if ( OrderofHS == 2 )
                    Fnorm = norm ( ProjectionOnBallofBsq( : , : , n ), 'fro' ) ;
                    if ( Fnorm > 1 )
                        ProjectionOnBallofBsq( : , : , n ) = ProjectionOnBallofBsq( : , : , n ) / Fnorm ;
                    end
                
                % q = inf , OrderofHS = p =1, 1/p + 1/q = 1
                %     [ Usvd, ProjectionOnBallofBsq, Vsvd ] = svd ( ProjectionOnBallofBsq( : , : , n ) ) ; 
                end
                
          end %n pixel
          Knesterovcurrent = ( 1 + sqrt( 1 + Knesterovprevious^2 ) ) / 2 ;
          Omegacurrent = ProjectionOnBallofBsqprevious + ( Knesterovprevious - 1 ) / Knesterovcurrent .* ( ProjectionOnBallofBsq - ProjectionOnBallofBsqprevious ) ;
          % stop condition
          if ( StopDeterminer( ProjectionOnBallofBsq , ProjectionOnBallofBsqprevious , inner_threshold ) )
              break ;
          end
          Knesterovprevious = Knesterovcurrent ;
          Omegaprevious = Omegacurrent;         
    end  %nesterov
    
    S = Z - Tao * HstarOmega( ProjectionOnBallofBsqprevious , ps) ;
    S_single = single( S ) ;
    AS = ProjectionCone_3D (S_single, t_length, s_length, z_length, Size, BetaScanRange, Pdomain, Xigamadomain, Distance) ;
    Kmmcurrent = ( 1 + sqrt( 1 + Kmmprevious^2 ) ) / 2 ;
    
    % determine whether the obective function descents in this round of
    % iteration
    Smatrix = reshape( S , t_length, s_length, z_length ) ; 
    Cfaicurrent = ObjectiveFunction( V , AS , Smatrix, OrderofHS, Tao) ;
    if ( Cfaicurrent > Cfaiprevious )
        Ucurrent = Uprevious ;
    else
        Ucurrent = S ;
    end
    Xcurrent = Ucurrent + ( Kmmprevious / Kmmcurrent ) * ( S - Ucurrent ) + ( Kmmprevious - 1) / Kmmcurrent * ( Ucurrent - Uprevious ) ;
    % stop condition
    if ( StopDeterminer( Ucurrent , Uprevious , outer_threshold ) )
              break ;
    end
    
    Uprevious = Ucurrent ; 
    Cfaiprevious = Cfaicurrent ;
    Xprevious = Xcurrent ; 
    Kmmprevious = Kmmcurrent ; 
end  %mm

Display = Ucurrent ; 
Display = reshape( Display , t_length, s_length, z_length ) ;
figure,imshow3Dfull ( Display , [] )

toc
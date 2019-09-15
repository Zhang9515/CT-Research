% 2017/06/25 by ZXZ
% modified on 2018/03/23
% 2018/04/09 apply GPU acceleration
% implement axial UIData
tic 
clear ;
Size = [ 512 ; 512 ; 23 ] ;     % actual range ( ScanFoV: 500 ) , 50 slice reconstruction images ( mm )
PicSize = 512 ; % that menas the resolution is 1 in the horizontal plane, and the vertical resolution is Resolution_z
t_length = PicSize ; s_length = PicSize ; z_length = 50 ;              % store the size of picture, z_length is 2 * 39.5*0.915/1060.2*( 570 - 256 ) = 23
Resolution_z = Size(3) / z_length ;
Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = Size(3) / 2 ;          % define the center 

Distance = 570 ;                                                                                              % distance between source and center point ( mm )
Distance_s2d = 1060.2 ;                                                                                       % distance between source and detector ( mm )

LGama =  936 ;                                                                                     
GamaInt = 9.7527e-4 ;          % Gama ( radian ) 
fid = fopen( 'G:\CTcode\Data\UIData\1\detectorangle.dat' , 'rb' ) ;         
[ Gamadomain ,  ~ ] = fread ( fid , LGama , 'float' ) ;
fclose ( fid ) ;
% GamaMax = 0.4526 ;                                                          % half max fan angle
GamaMax = abs ( Gamadomain ( end ) - Gamadomain ( 1 ) ) ;
Gamadomain = -Gamadomain ;

XigamaInt = 0.9150 ;                                      % interval of Xigama ( 0.5 acq-SliceThickness or 0.5508 acq-ChanSpace or DetWidth 0.9150 )
MaxXigama = 39.5 * XigamaInt ;
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = Xigamadomain' ;
LXigama = length ( Xigamadomain ) ;
 
BetaScanInt = deg2rad ( 240/800 ) ;             % scanning internal   ( Angle )  ScanArc: 240
MaxBeta = deg2rad ( 240 ) ;
StartAngle =  deg2rad ( 1.5466 ) ; 
BetaScanRange = StartAngle : BetaScanInt : ( MaxBeta + StartAngle - BetaScanInt ) ;     % scanning range , angle between SO and aixs Y
BetaScanRange = BetaScanRange' ;
LBeta = length ( BetaScanRange ) ; 

fid = fopen ( 'G:\CTcode\Data\UIData\1\data.dat' , 'rb' ) ;                                   % load projection data
[ R , count ] = fread ( fid , LGama * LXigama * 800 , 'float') ;                   % UI data1 all have value, but data2 only beta 437-800 have value
% R = reshape ( R , LGama , LXigama , 800 ) ;
%R = R ( LGama * LXigama * 436 + 1: end );
fclose ( fid ) ;

%  this is wrong, actually it should inverse the sign of Gama
% inverse the order of Gama ********* important ***********
% for Betaindex = 1 : LBeta
%         for Xigamaindex  = 1 : LXigama
%                 Media = R ( : , Xigamaindex , Betaindex ) ;
%                 R ( : , Xigamaindex , Betaindex ) = flipud ( Media ) ;
%         end
% end
%% pre-weighting and forming SL
% G = zeros ( 2 * LGama - 1  , LXigama ) ;                                                      % create filter
% R = reshape ( R , LGama , LXigama , 800 ) ;
% for Gamaindex = 1 : LGama
%     for Xigamaindex  = 1 : LXigama
%         
%         Proportion = Distance2 / sqrt ( Distance2^2 + Xigamadomain( Xigamaindex )^2 ) ;
%         gama_1 =  Gamadomain ( Gamaindex ) * Proportion ;       
%         
%         gama =  Gamadomain ( Gamaindex ) ;      
%         
%         Distance_1 = Distance / Proportion ; 
%         for Betaindex = 1 : LBeta
%             
%             beta = BetaScanRange ( Betaindex ) ;
%             % short scan weight function
% %                 beta_1 = BetaScanRange ( Betaindex ) * Proportion ;
% %                     GamaMax_1 = GamaMax * Proportion ;
% %             if ( beta >= StartAngleRadian && beta < StartAngleRadian + 2 * GamaMax - 2 * gama )
% %                 W = sin ( pi / 4 * beta_1 / ( GamaMax_1 - gama_1 ) )^2 ;
% %             elseif ( beta >= StartAngleRadian + 2 * GamaMax - 2 * gama && beta <= StartAngleRadian + pi - 2 * gama )
% %                 W =  1 ;
% %             elseif ( beta > StartAngleRadian + pi - 2 * gama && beta <= StartAngleRadian + pi + 2 * GamaMax )
% %                 W = sin ( pi / 4 * ( pi + 2 * GamaMax_1 - beta_1 ) / ( GamaMax_1 + gama_1 ) )^2 ;
% %             else 
% %                 W = 0 ;
% %             end
%           
%             W = ParkFunction ( beta , StartAngle , gama , GamaMax ) ;
%             
%             R ( Gamaindex , Xigamaindex , Betaindex ) = R ( Gamaindex , Xigamaindex , Betaindex ) * Distance * cos ( gama ) * W ; 
%             
%         end
%         %%% form S_L filter %%%%%%%%
%         %%%%%%%%%%%%%%%%%%%
% %             Gamaindex_SL =  Gamaindex - 1 ;
% %             if Gamaindex_SL == 0
% %                 G ( LGama , Xigamaindex ) = 2 / ( pi * GamaInt )^2 ;    % radian 
% %             else                                       
% %                G ( LGama + Gamaindex_SL , Xigamaindex ) = - 2 * Gamaindex_SL^2 / pi^2 / ( 4 * Gamaindex_SL^2 - 1 ) / ( sin ( Gamaindex_SL * GamaInt ) )^2 ;
% %                G ( LGama - Gamaindex_SL , Xigamaindex ) = - 2 *  Gamaindex_SL^2 / pi^2 / ( 4 * Gamaindex_SL^2 - 1 ) / ( sin ( Gamaindex_SL * GamaInt ) )^2  ;
% %             end
%         
%     end
% end
% % figure , imshow ( squeeze ( R ( : , 40  , : ) ) , [] ) ; 
% %  save R0401.mat R 
% 
%% filter    S-L
% 
% Rcov = zeros ( LGama , LXigama , LBeta ) ;
% parfor Betaindex = 1 :  LBeta
%     for Xigamaindex = 1 : LXigama
%         
%         Proportion = Distance2 / sqrt ( Distance2^2 + Xigamadomain ( Xigamaindex )^2 ) ;
%         GamaInt_1 = GamaInt * Proportion ;
%         cov = conv ( squeeze ( R ( : , Xigamaindex , Betaindex ) ) , G ( : , Xigamaindex )' ) ;                              % convolution with filter
%         Rcov ( : , Xigamaindex , Betaindex ) = GamaInt_1 * cov ( LGama : 2 * LGama - 1 ) ;
% 
%     end
% end
% 
% Rcov = zeros ( LGama , LXigama , LBeta ) ;
% parfor Betaindex = 1 :  LBeta
%     for Xigamaindex = 1 : LXigama
%         
%         %Proportion = Distance2 / sqrt ( Distance2^2 + Xigamadomain ( Xigamaindex )^2 ) ;
%         %GamaInt_1 = GamaInt * Proportion ;
%         cov = conv ( squeeze ( R ( : , Xigamaindex , Betaindex ) ) , G ( : , Xigamaindex )' ) ;                              % convolution with filter
%         Rcov ( : , Xigamaindex , Betaindex ) = GamaInt * cov ( LGama : 2 * LGama - 1 ) ;
% 
%     end
% end
% % save R04011.mat R 
% toc
% % % figure , imshow ( squeeze ( Rcov ( : , 40  , : ) ) , [] ) ; 
% % clear R G ;
%% reconstruct
% 
% Display = zeros ( PicSize , PicSize ,  z_length ) ;
% ProjectionSlice = zeros ( PicSize , PicSize , LBeta ) ;
% ProjectionSlice = reshape ( ProjectionSlice , 1 , [] ) ;
% for z = 1 : z_length
% 
% %      z = 1 ;
%     
%     for num = 1 : t_length * s_length * LBeta
% 
%                t = floor ( ( num - 1 ) / ( LBeta * s_length ) ) + 1 ;
% 
%                s = floor ( ( num - ( t - 1 ) * s_length * LBeta - 1 ) / LBeta ) + 1 ; 
% 
%                betas = mod ( num - 1 , LBeta ) + 1 ;
%                
%                     beta =  BetaScanRange ( betas ) ;
%                     
%                     source_t = Distance * cos ( beta + pi / 2 ) ;                           % define the source in ground coordinate
%                     source_s = Distance * sin ( beta + pi / 2 ) ;   
%                     source_z = 0 ; 
%                                         
%                     image_t = ( t - 0.5 ) - Center_t  ;  image_s = ( s - 0.5 ) - Center_s  ; image_z = ( z - 0.5 ) * Resolution_z - Center_z  ;           % image pixel in ground coordinate
%                     L2 = ( image_t - source_t )^2 + ( image_s - source_s )^2 + ( image_z - source_z )^2 ;
%                     
%                     dect_t = image_t * cos ( beta ) + image_s * sin ( beta ) ;          % rotate in ground coordinate
%                     dect_s = - image_t * sin ( beta ) + image_s * cos ( beta ) ; 
%                     dect_z = image_z ;     
%                     
%                     Xigama = Distance2 * tan ( asin ( dect_z / sqrt ( L2 ) ) ) ;                 % define the projection position on the detector       
%                     Gama = atan ( dect_t / ( Distance - dect_s ) ) ; 
%                     
%                     Proportion = sqrt ( L2 ) / sqrt ( L2 - dect_z^2 ) ;
%                     if  ( Gama >= min ( Gamadomain ) && Gama < max ( Gamadomain )  && Xigama >= min ( Xigamadomain ) && Xigama < max ( Xigamadomain ) )
%                             
%                         XigamaN1index =  floor ( abs ( Xigama - Xigamadomain( 1 ) ) / XigamaInt ) + 1 ; 
%                         XigamaN2index = XigamaN1index + 1 ;
%                         GamaN1index =  floor ( abs ( Gama - Gamadomain( 1 ) ) / GamaInt ) + 1 ; 
%                         GamaN2index = GamaN1index + 1 ;
%                         
%                             Gama_domain1 = Gamadomain ( GamaN1index ) ; Gama_domain2 = Gamadomain ( GamaN2index ) ;
%                             Xigama_domain1 = Xigamadomain ( XigamaN1index ) ; Xigama_domain2 = Xigamadomain ( XigamaN2index ) ;  
%                             
%                             % bilinear interpolation
%                             Xig1 = abs ( Xigama - Xigama_domain1 ) ; Xig2 = abs ( Xigama_domain2 - Xigama )  ;
%                             Gama1 = abs ( Gama - Gama_domain1 ) ; Gama2 = abs ( Gama_domain2 - Gama ) ;
%                             ProjectionSlice ( num ) = ( Xig2 * Gama2 * Rcov ( GamaN1index , XigamaN1index , betas ) ...
%                                 + Xig1 * Gama2 * Rcov ( GamaN1index , XigamaN2index , betas ) + Xig2 * Gama1 * Rcov ( GamaN2index , XigamaN1index , betas ) ...
%                                 + Xig1 * Gama1 * Rcov ( GamaN2index , XigamaN2index , betas ) ) / ( GamaInt * XigamaInt ) / L2 *  BetaScanInt ;
%                     end                       
%     end
%     
%      for  t = 1 : t_length 
%             for s = 1 : s_length 
%                     Display ( t , s , z ) = sum ( ProjectionSlice ( 1 + ( ( t - 1 ) * s_length * LBeta + ( s - 1 ) * LBeta ) : LBeta + ( ( t - 1 ) * s_length * LBeta + ( s - 1 ) * LBeta ) ) ) ;
%             end
%      end   
% 
% end

%% ImplementUI CUDA GPU accelerate

Display = ImplementUI ( R , Xigamadomain , Gamadomain , BetaScanRange ,  StartAngle , GamaMax , Distance , Distance_s2d , Size , z_length ) ;
Display = reshape ( Display , t_length , s_length , z_length ) ;

% clear Rcov ProjectionSlice ;
% figure , imshow ( flipud ( Display ) , [ 100,2000 ] ) ; 

toc ;
%% display

% a = open ( 'D:\TestCpp\CT\Data\UIData\1\Rec24.fig' ) ;
% h = get (gca , 'Children') ;
% data = get ( h , 'Cdata') ;
% figure , imshow ( data, [ 100,2000] ) ; 







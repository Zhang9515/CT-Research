% 2019/01/10 by ZXZ
% GPUFDK, all codes are extracted from threeDprojection.m, just to make the
% file clear and simple
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
clear;

load('E:\ZXZ\Data\ThoraxHD_661.mat')
% load('G:\CTcode\Data\ThoraxHD_661.mat')
pic = ThoraxHD_661 ;
pic = single ( pic ) ;
clear ThoraxHD
% pic = single( Diskphantom ( 512 ) ) ;
% figure,imshow3Dfull( pic( : , : , 129 : 384 ) , [0 0.5] )

[ t_length , s_length , z_length ] = size ( pic ) ;              % store the size of picture ( pixels )

% Size = [ t_length * 0.7480 ; s_length * 0.7480 ; z_length * 1 ] ;     % actual range ( should refer to dicom info)
Size = [ 60 , 60 , 60 * z_length / t_length ] ;

z_lengthrec = z_length ;    % backprojection matrix size ( usually same as the size of original image) 
Sizerec = [ 60 , 60 , 60 * z_lengthrec / t_length ] ;     % reconstruction real size

Rplane = Size(1) * sqrt ( 2 ) / 2 ;                    %  circumscribed circle radius of project in the plane

Resolution = [ Size(1) / t_length , Size(2) / s_length , Size(3) / z_length ];                           

Distance = 207.6 ;              % distance between source and center point ( 207.6 for 60 )

PInt = 0.2 ;                                    %interval of P ( 0.1 exact )
MaxP = max ( ( Size(1) / 2 * Distance ) / ( Distance - Size(1) / 2 ) , Rplane * 1.1 );
Pdomain = - MaxP : PInt : MaxP ;                       % detective range P
Pdomain = single(Pdomain');
LP = length ( Pdomain ) ;

XigamaInt = 0.2 ;                                      % interval of Xigama ( 0.1 exact )
MaxXigama = ( Sizerec(3) / 2 * Distance ) / ( Distance - Size(1) * sqrt(2) / 2 ) ;       % computed by rule of similar triangle
MaxXigama = 46.669 ; 
Xigamadomain = - MaxXigama : XigamaInt : MaxXigama ;                       % detective range Xigama
Xigamadomain = single(Xigamadomain');
LXigama = length ( Xigamadomain ) ;

Center_t = Size(1) / 2 ;  Center_s = Size(2) / 2 ;   Center_z = Size(3) / 2 ;          % define the center 

% FanAmax = atan ( MaxP / ( Distance - Rplane ) ) ; 
BetaScanInt = deg2rad(4) ;             % scanning internal    ( 0.3 exact )           
% MaxBeta = 180 + 2 * FanAmax * 180 / pi ;         % short scan 
MaxBeta = deg2rad(360) ;
BetaScanRange = BetaScanInt : BetaScanInt : MaxBeta  ;     % scanning range , angle between SO and aixs Y
BetaScanRange = single(BetaScanRange');
LBeta = length ( BetaScanRange ) ; 

% to ensure the projection matrix not too big, here compute LP*LXigama*LBeta to
% control its value less than 400 * 400 * 90. if not, split the Beta range
% separately 

picvector = reshape (pic, t_length * s_length * z_length, 1);
% clear pic;

Display = zeros ( t_length * s_length * z_lengthrec, 1 ) ; 
LBetaPrime = floor ( ( 700 * 700 * 90 ) / ( LP * LXigama ) ) ;     % 700 * 700 * 90 is related to the quality of GPU(1080ti)
                                                                                                   % 300*300*90 (1060)
times = ceil ( LBeta / LBetaPrime ) ;
disp(['times: ',num2str(times)]);

for  split = 1 : times 
    
    disp(['split: ',num2str(split)]);
    if ( split ~= times )       
        BetaScanRangePrime = BetaScanRange ( 1 + ( split -1 ) * LBetaPrime : split * LBetaPrime ) ;
    else
        LBetaPrime = LBeta - (split  - 1 ) * LBetaPrime ;
        BetaScanRangePrime = BetaScanRange ( 1 + ( split -1 ) * LBetaPrime : LBeta ) ;
    end
    
    %% GPU accelerate projection

    % size of R is ( LP * LXigama * LBeta, 1 )
    R = ProjectionCone_3D (picvector, t_length, s_length, z_length, Size, BetaScanRangePrime, Pdomain, Xigamadomain, Distance);
%     R = reshape( R , LP , LXigama , LBetaPrime ) ;
%     figure,imshow3Dfull(R, [], 'grey')

    %% GPU accelerate the preweight & filteration & backprojection

     R= single(R) ;
    Display = Display + Backprojection ( R , Xigamadomain , Pdomain , BetaScanRangePrime , Distance, Sizerec, t_length, s_length, z_lengthrec) ;
    clear R ;
end
Display = reshape ( Display , t_length , s_length , z_lengthrec ) ;
clear pic picvector
% figure,imshow3Dfull(Display, [], 'grey')

toc


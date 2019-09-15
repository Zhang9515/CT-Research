% 2016/12/20
% ZXZ
% according to 'Studies on implementation of the Katsevich' by Hengyong Yu and Ge Wang
%
%                                /
%                              / 
%                            /
%     -----------------/------------------>  x2
%                       /
%                    /
%               | /     x1
%                                                    x1 * x2 = x3

clear all ;
SScanInt = 1 ;             % scanning internal              
MaxS = 180 * 3 ; 
SScanRange =  -MaxS : SScanInt : MaxS ;     % scanning range , angle between SO and aixs S
pic = phantom3d ( 32 ) ;     % original picture  
[ X1_length , X2_length , X3_length ] = size ( pic ) ;              % store the size of picture
r = max ( size ( pic ) ) * sqrt ( 2 ) / 2 ;                                      % radius of object
RRatio = 2 ; 
Radius = max ( size ( pic ) ) * RRatio ;                                     % distance between source and center point
miu0 = r / Radius ;                                                                    % ratio  between radius of object and scanning trajectory ( nearly 1/3 )
DeltaC = 2 * acos ( r / Radius ) ;                                                       % value of Delta
DisRatio = 2 * RRatio ;  
Distance = max ( size ( pic ) ) * DisRatio ;                               % distance between source and detect plain
HelicalPitch = max ( size ( pic ) ) * 0.5 ;                                  % height per turn of helix
UInt = 1 ;                                                                                   %interval of U
MaxU = round ( Distance * sin ( DeltaC ) * ( 1 - cos ( DeltaC ) ) ) ;
Udomain = - MaxU : UInt : MaxU ;                                         % detective range U
VInt = 1 ;                                                                                   %interval of V
MaxV = round ( Distance * HelicalPitch * ( 2 * pi -  DeltaC ) / ( 2 * pi * Radius * ( 1 - cos ( DeltaC ) ) ) ) ;
Vdomain = - MaxV : VInt : MaxV ;                                          % detective range V
Center_X1 = X1_length / 2 ;  Center_X2 = X2_length / 2 ;   Center_X3 = X3_length / 2 ;          % define the center 

%% projection
% define position of center point in the helix : ( R , 0 , 0 )
R = zeros ( length ( SScanRange ) ,  length ( Udomain ) , length ( Vdomain ) ) ;   % create space to store cone projection
distanceInt = 1 ;      %define the interval on the rayline
for i = 1 : length ( SScanRange )
    s =  SScanRange ( i ) ;
    SRadian = s * pi / 180 ;                    % angle to radian
    source_X1 = Center_X1 + Radius * cos ( SRadian ) ; source_X2 = Center_X2 + Radius * sin ( SRadian ) ; source_X3 = Center_X3 + SRadian * HelicalPitch / ( 2 * pi ) ; 
    for j = 1 : length ( Vdomain )           
        theta = atan ( Vdomain ( j ) / Distance ) ;          % radian angle in X1'-X2 coordinate plane
        Distance_shift = Distance / cos ( theta ) ;                    % length of DO'
        for k = 1 : length ( Udomain )
                gama = atan ( Udomain ( k ) / Distance_shift ) ;                                    % radian angle in X1''-X2 coordinate plane 
           
            %   Siddon algrithim
                Smax = 2 * Distance ; 
                DetectPoint_X1end = Center_X1 + ( Radius  - Smax * cos ( theta ) * cos ( gama ) ) * cos ( SRadian ) - Smax * sin ( SRadian ) * sin ( gama ) ; 
                DetectPoint_X2end = Center_X2 + Smax * sin ( gama ) * cos ( SRadian ) + ( Radius  - Smax * cos ( theta ) * cos ( gama ) ) * sin ( SRadian ) ;
                DetectPoint_X3end = Center_X3 + SRadian * HelicalPitch / ( 2 * pi ) + Smax * cos ( gama ) * sin ( theta ) ;
                if ( DetectPoint_X1end > source_X1 )                                                         % define projection line direction and range of i and j
                     X1_signal = 1 ;   
                     X1_range = ceil ( source_X1 ) : floor ( DetectPoint_X1end ) ;
                else 
                    X1_signal = -1 ; 
                    X1_range = ceil ( DetectPoint_X1end ) : floor ( source_X1 ) ;
                end
                if ( DetectPoint_X2end > source_X2 )
                    X2_signal = 1 ; 
                    X2_range = ceil ( source_X2 ) : floor ( DetectPoint_X2end ) ;
                else
                    X2_signal = -1 ; 
                    X2_range = ceil ( DetectPoint_X2end ) : floor ( source_X2 ) ;
                end     
                if ( DetectPoint_X3end > source_X3 )
                    X3_signal = 1 ; 
                    X3_range = ceil ( source_X3 ) : floor ( DetectPoint_X3end ) ;
                else
                    X3_signal = -1 ; 
                    X3_range = ceil ( DetectPoint_X3end ) : floor ( source_X3 ) ;
                end
                LX1 = length ( X1_range ) ; LX2 = length ( X2_range ) ; LX3 = length ( X3_range ) ;
                ProjectionLine = zeros ( 3 , LX1 + LX2 + LX3 ) ;
                ProjectionLine ( 1 , 1 : LX1 ) = abs ( ( X1_range - source_X1 ) / ( DetectPoint_X1end - source_X1 + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , 1 : LX1 ) = X1_range ;
                ProjectionLine ( 3 , 1 : LX1 ) = 1 ;                              % 1 represents t
                ProjectionLine ( 1 , ( LX1 + 1 ) : ( LX1 + LX2 ) ) = abs ( ( X2_range - source_X2 ) / ( DetectPoint_X2end - source_X2 + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( LX1 + 1 ) : ( LX1 + LX2 ) ) = X2_range ;
                ProjectionLine ( 3 , ( LX1 + 1 ) : ( LX1 + LX2 ) ) = 2 ;                              % 2 represents s
                ProjectionLine ( 1 , ( LX1 + LX2 + 1 ) : ( LX1 + LX2 + LX3 ) ) = abs ( ( X3_range - source_X3 ) / ( DetectPoint_X3end - source_X3 + 1e-10 ) ) * Smax ;
                ProjectionLine ( 2 , ( LX1 + LX2 + 1 ) : ( LX1 + LX2 + LX3 ) ) = X3_range ;
                ProjectionLine ( 3 , ( LX1 + LX2 + 1 ) : ( LX1 + LX2 + LX3 ) ) = 3 ;                              % 3 represents z
                ProjectionLine = sortrows ( ProjectionLine' )' ;
                DetectPoint_X1 = ceil ( source_X1 ) ; DetectPoint_X2 = ceil ( source_X2 ) ; DetectPoint_X3 = ceil ( source_X3 ) ;                   % start point of the line        
                
                for n = 1 : length ( ProjectionLine )
                    
                    if ( ProjectionLine ( 3 , n ) == 1 )
                            DetectPoint_X1 = ProjectionLine ( 2 , n ) + 0.5 +  0.5 * X1_signal ;                    % up/right is 1 , down/left is 0
                    elseif ( ProjectionLine ( 3 , n ) == 2 )
                            DetectPoint_X2 = ProjectionLine ( 2 , n ) + 0.5 +  0.5 * X2_signal ; 
                    elseif ( ProjectionLine ( 3 , n ) == 3 )
                            DetectPoint_X3 = ProjectionLine ( 2 , n ) + 0.5 +  0.5 * X3_signal ; 
                    end
                    
                    if ( LX1 == 0 && source_X1 == 0 && DetectPoint_X2 > 0 && DetectPoint_X2 <= X2_length ...                          % boundary condition
                            && DetectPoint_X3 > 0 && DetectPoint_X3 <= X3_length )   
                        
                        f = pic ( 1 , DetectPoint_X2 , DetectPoint_X3 ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , k , j ) = R ( i , k , j ) + f * LengthUnit ;    
                    end    
                    if ( LX2 == 0 && source_X2 == 0 && DetectPoint_X1 > 0 && DetectPoint_X1 <= X1_length ...
                            && DetectPoint_X3 > 0 && DetectPoint_X3 <= X3_length )   
                        
                        f = pic ( DetectPoint_t , 1 , DetectPoint_z ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , k , j ) = R ( i , k , j ) + f * LengthUnit ;    
                    end      
                    if ( LX3 == 0 && source_X3 == 0 && DetectPoint_X1 > 0 && DetectPoint_X1 <= X1_length ...
                            && DetectPoint_X2 > 0 && DetectPoint_X2 <= X2_length )    
                        
                        f = pic ( DetectPoint_X1 , DetectPoint_X2 , 1 ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( i , k , j ) = R ( i , k , j ) + f * LengthUnit ;    
                    end   
                    
                    if ( DetectPoint_X1 > 0 && DetectPoint_X1 <= X1_length && DetectPoint_X2 > 0 && DetectPoint_X2 <= X2_length && DetectPoint_X3 > 0 && DetectPoint_X3 <= X3_length ) 
                         
                        if ( LX1 == 0 && source_X1 == DetectPoint_X1 )
                             if ( source_X1 == X1_length )
                                 f = pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) + pic ( DetectPoint_X1 + 1 , DetectPoint_X2 , DetectPoint_X3 ) ) / 2 ;
                             end
                         elseif ( LX2 == 0 && source_X2 == DetectPoint_X2 ) 
                             if ( source_X2 == X2_length )
                                 f = pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) + pic ( DetectPoint_X1 , DetectPoint_X2 + 1 , DetectPoint_X3 ) ) / 2 ;
                             end
                         elseif ( LX3 == 0 && source_X3 == DetectPoint_X3 ) 
                             if ( source_X3 == X3_length )
                                 f = pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) + pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 + 1 ) ) / 2 ;
                             end
                         else 
                                f = pic ( DetectPoint_X1 , DetectPoint_X2 , DetectPoint_X3 ) ;
                         end
                  
                         LengthUnit =  ( ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ) ; 
                         R ( i , k , j ) = R ( i , k , j ) + f * LengthUnit ;    
                    end
                end 
                       
        end
    end 
end
% figure , imshow ( squeeze ( R ( 361 , : , : ) ) , [] ) ; 
% 
% save D:\TestCpp\R.mat R ;
% load D:\TestCpp\R.mat ;

%% derivation 
RDeri = zeros ( length ( SScanRange ) ,  length ( Udomain ) , length ( Vdomain ) ) ;   % create space to store fan projection
for i = 1 : length ( SScanRange )
    for j = 1 : length ( Udomain )
        for k = 1 : length ( Vdomain )
            
            if ( i == 1 )
                RsDeri = ( R ( i + 1 , j , k ) - R ( i , j , k ) ) / ( SScanInt * pi / 180 ) ;
            elseif ( i == length ( SScanRange ) )
                RsDeri = ( R ( i , j , k ) - R ( i - 1 , j , k ) ) / ( SScanInt * pi / 180 ) ;
            else 
                RsDeri = ( R ( i + 1 , j , k ) - R ( i - 1 , j , k ) ) / ( 2 * SScanInt * pi / 180 ) ;
            end
            
            if ( j == 1 )
                RuDeri = ( R ( i , j + 1 , k ) - R ( i , j , k ) ) / UInt ;
            elseif ( j == length ( Udomain ) )
                RuDeri = ( R ( i , j , k ) - R ( i , j - 1 , k ) ) / UInt ;
            else 
                RuDeri = ( R ( i , j + 1 , k ) - R ( i , j - 1 , k ) ) / ( 2 * UInt ) ;
            end
            
            if ( k == 1 )
                RvDeri = ( R ( i , j , k + 1 ) - R ( i , j , k ) ) / VInt ;
            elseif ( k == length ( Vdomain ) )
                RvDeri = ( R ( i , j , k ) - R ( i , j , k - 1 ) ) / VInt ;
            else 
                RvDeri = ( R ( i , j , k + 1 ) - R ( i , j , k - 1 ) ) / ( 2 * VInt ) ;
            end
            RDeri ( i , j , k ) = RsDeri + RuDeri * ( Distance ^ 2 + Udomain ( j ) ^ 2 ) / Distance + RvDeri * Udomain ( j ) * Vdomain ( k ) / Distance ; 
        end
     end
end
% save D:\TestCpp\RDeri.mat RDeri ;
% load  D:\TestCpp\RDeri.mat ;

%% filtration 
% r = R , hence, range of Sq : [ - pi + 0.5dt , pi - 0.5dt ]
SqScanInt = 3 ;             % scanning internal              
MaxSq = ( pi - 0.5 * DeltaC ) * 180 / pi   ;            % S1max
MinSq = ( -pi + 0.5 * DeltaC ) * 180 / pi  ;            % S1min
SqScanRange =  MinSq : SqScanInt : MaxSq;     % scanning range to form the filtration plane

% compute Vq
Vqdomain = zeros ( length ( SqScanRange ) , length ( Udomain ) ) ;   % store the v value in k-line 
for i = 1 : length ( SqScanRange )
    for j = 1 : length ( Udomain )
        SQ = SqScanRange ( i ) * pi / 180 ;
        if ( SqScanRange ( i ) ~= 0 )
            VV = HelicalPitch * SQ / ( 2 * pi * Radius ) * ( Distance + cot ( SQ ) * Udomain ( j ) ) ; 
        else 
            VV = HelicalPitch * Udomain ( j ) / ( 2 * pi * Radius ) ; 
        end
        Vqdomain ( i , j ) = VV ;
    end
end
% save D:\TestCpp\Vqdomain.mat Vqdomain ;
% load D:\TestCpp\Vqdomain.mat ;

% compute fai_1
fai1 = zeros ( length ( SScanRange ) , length ( SqScanRange ) , length ( Udomain ) ) ;              % create space to store fan projection
for i = 1 : length ( SScanRange )
    for j = 1 : length ( SqScanRange )
        for k = 1 : length ( Udomain )
                [ ~ , Vqmindex ] = min ( abs ( Vdomain - Vqdomain ( j , k ) ) ) ;
                fai1 ( i , j , k ) = RDeri ( i , k , Vqmindex ) * Distance / sqrt ( Distance ^ 2 + Udomain ( k ) ^ 2 + Vqdomain ( j , k ) ^ 2 ) ;
        end
    end
end

% save D:\TestCpp\fai1.mat fai1 ;
% load D:\TestCpp\fai1.mat ;

% compute fai_2 and fai_3
fai23 = zeros ( length ( SScanRange ) ,  length ( SqScanRange ) , length ( Udomain ) ) ;   % create space to store fan projection
for i = 1 : length ( SScanRange )
    for j = 1 : length ( SqScanRange )
        for k = 1 : length ( Udomain )
            for h = 1 : length ( Udomain )
                if ( h ~= k )
                    fai23 ( i , j , k ) = fai23 ( i , j , k ) + UInt * fai1 ( i , j , h ) / ( Udomain ( h ) - Udomain ( k ) ) ; 
                end             
            end
            fai23 ( i , j , k ) = fai23 ( i , j , k ) * sqrt ( Distance ^ 2 + Udomain ( k ) ^ 2 + Vqdomain ( j , k ) ^ 2 ) / Distance ;
        end
    end
end

% save D:\TestCpp\fai23.mat fai23 ;
% load D:\TestCpp\fai23.mat ;

% compute fai interpolated by fai3
fai = zeros ( length ( SScanRange ) ,  length ( Udomain ) , length ( Vdomain ) ) ;   % create space to store fan projection
for i = 1 : length ( SScanRange )
    for j = 1 : length ( Udomain )
        for k = 1 : length ( Vdomain )
            [ ~ , Vqmindex ] = min ( abs ( Vqdomain ( : , j ) - Vdomain ( k ) ) ) ;
            MinVq = min ( Vqdomain ( : , j ) ) ;
            MaxVq = max ( Vqdomain ( : , j ) ) ;
            if ( Vdomain ( k ) >= MinVq  &&  Vdomain ( k ) <= MaxVq )
                fai ( i , j , k ) = fai23 (  i , Vqmindex , j ) ;
            else
                fai ( i , j , k ) = 0 ;
            end
        end
    end
end

% save D:\TestCpp\fai.mat fai ;

% load D:\TestCpp\fai.mat ;

%% backprojection

Projection = zeros ( X1_length , X2_length , X3_length ) ;      % reconstruct the pic
for i = 1 : X1_length
    for j = 1 : X2_length
        for k =1 : X3_length
            % compute two peak points of pi-line 
            X10 = i - Center_X1 - 0.5 ; X20 = j - Center_X2 - 0.5 ; X30 = k - Center_X3 - 0.5 ;            
            S0 = 2 * pi * X30 / HelicalPitch ; 
            if ( X20 > 0 )
                xigama =  atan ( X20 / ( X10 + 1e-10 )  ) ;   
            else
                xigama = -atan ( X20 / ( X10 + 1e-10 ) ) ;
            end
            r0 = sqrt ( X10^2 + X20^2 ) ;
            % Sb search dom ain , compute by iteration in one-variable equation
            Sb0max = S0 - acos ( miu0 ) * ( 1 - miu0 );
            Sb0min = S0 - ( pi - acos ( miu0 ) ) * ( 1 + miu0 ) ;
            Sb0 = Sb0min + ( Sb0max - Sb0min ) / 2 ; 
               % iteration 
%             iteration_count = 0 ;                                                % iteration times 
%             while (1)
%                 iteration_count = iteration_count + 1 ;
%                 t = ( Radius^2 + r^2 ) / ( 2 * Radius * ( Radius - r * cos ( xigama - Sb0 ) ) ) ;
%                 Sb1 = S0 + 2 * ( t - 1 ) * acos ( r * sin ( xigama - Sb0 ) / sqrt ( Radius^2 + r^2 - 2 * Radius * r * cos ( xigama - Sb0 ) ) ) ; 
%                 err = abs ( ( Sb1 - Sb0 ) * 180 / pi ) ;                   %  iteration condition
%                 if ( err  < 0.5  )
%                     break ; 
%                 elseif ( iteration_count > 100 )
%                     disp(' cant solve Sb ') ;
%                     return ;
%                 end 
%                 Sb0 = Sb1 ;
%             end          

            % dichotomy
            Sb0maxAng = Sb0max * 180 / pi ;
            Sb0minAng = Sb0min * 180 / pi ; 
            while ( 1 )
                if ( Sb0maxAng - Sb0minAng < 0.5 )
                    Sb1 = Sb0 ;
                    break ;
                end
                Sb0Ang = Sb0minAng + ( Sb0maxAng - Sb0minAng ) / 2 ; 
                Sb0 = Sb0Ang * pi / 180 ;
                t = ( Radius^2 - r0^2 ) / ( 2 * Radius * ( Radius - r0 * cos ( xigama - Sb0 ) ) ) ;
                St = Sb0 + 2 * acos ( r0 * sin ( xigama - Sb0 ) / sqrt ( Radius^2 + r0^2 - 2 * Radius * r0 * cos ( xigama - Sb0 ) ) ) ;
                X31 = HelicalPitch / ( 2 * pi ) * ( Sb0 * t + St * ( 1 - t ) ) ;
                if ( X31 > X30 )
                    Sb0maxAng = Sb0Ang ; 
                elseif ( X31 < X30 )
                    Sb0minAng = Sb0Ang ;
                else 
                    Sb1 = Sb0 ;
                    break ;
                end
            end
            St = round ( ( Sb1 + 2 * acos ( r0 * sin ( xigama - Sb1 ) / sqrt ( Radius^2 + r0^2 - 2 * Radius * r0 * cos ( xigama - Sb1 ) ) ) ) * 180 / pi ) ;    % obtain St
            Sb = round ( Sb1 * 180 / pi ) ;                                  % obtain Sb ( in angle domain )
            Sbindex = find ( SScanRange == Sb ) ; 
            Stindex = find ( SScanRange == St ) ; 
            % backprojection
            for Sindex = Sbindex : Stindex   
                Sradian = SScanRange ( Sindex ) * pi / 180 ;
                CosS = cos ( Sradian ) ;
                SinS  = sin ( Sradian ) ;
                d1 = [ -SinS , CosS ,  0 ] ;
                d2 = [ 0 , 0 , 1 ] ;
                d3 = [ -CosS , -SinS , 0 ] ;
                Ys3 = Sradian * HelicalPitch / ( 2 * pi ) ;
                Ys = [ Radius * CosS , Radius * SinS , Ys3 ] ; 
                delta = [ X10 , X20 , X30 ] - Ys ;
                UU = Distance * dot ( delta , d1 ) / dot ( delta , d3 ) ; 
                VV = Distance * dot ( delta , d2 ) / dot ( delta , d3 ) ; 
                UUindex1 = floor ( ( UU + MaxU ) / UInt ) + 1 ;
                UUindex2 = UUindex1 + 1 ; 
                VVindex1 = floor ( ( VV + MaxV ) / VInt ) + 1 ;
                VVindex2 = VVindex1 + 1 ; 
                if ( UU >= -MaxU && UU <= MaxU && VV >= -MaxV && VV <= MaxV )
                    UU1 = Udomain ( UUindex1 ) ; UU2 = Udomain ( UUindex2 ) ;
                    VV1 = Vdomain ( VVindex1 ) ; VV2 = Vdomain ( VVindex2 ) ;
                    
                    % Bilinear interpolation
                    UUd1 = UU - UU1 ; UUd2 = UU2 - UU ;
                    VVd1 = VV - VV1 ; VVd2 = VV2 - VV ;                   
                    Projection ( i , j , k ) = Projection ( i , j , k ) + ( UUd2 * VVd2 * fai ( Sindex , UUindex1 , VVindex1 ) ...
                        + UUd1 * VVd2 * fai ( Sindex , UUindex2 , VVindex1 ) + UUd2 * VVd1 * fai ( Sindex , UUindex1 , VVindex2 ) ...
                        + UUd1 * VVd1 * fai ( Sindex , UUindex2 , VVindex2 ) ) / norm ( delta ) / ( UInt * VInt ) ;
                end
            end
            Projection ( i , j , k ) = - ( SScanInt * pi / 180 ) / ( 2 * pi ^ 2 ) *  Projection ( i , j , k ) ;
        end
    end
end
figure , imshow( squeeze ( Projection ( : ,  : , Center_X3 ) ) ) ; 
% save D:\TestCpp\KatProjection.mat Projection ;






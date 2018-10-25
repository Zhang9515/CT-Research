% 2017/11/17 ZXZ
% 2018/04/09
% ART
% So slow & So inaccurate
% reshape set arrary first in column
tic
clear;
% close all;
%%
% parameter define
thetaint = deg2rad (0.5) ;                                                                 % theta unit 
thetaRange = thetaint : thetaint : pi ;                                % radon scanning range
thetaRange = thetaRange';                                               % GPU read column vector
Ltheta = length ( thetaRange ) ; 

pic = phantom ( 512 ) ;
Size = [ 60 , 60 ] ;                                  % actual range

[ height , width ] = size ( pic ) ;              % store the size of picture
Resolution = max ( Size ) / height ;   % define the resolution of pic
Center_y = Size ( 1 ) / 2 ;  Center_x = Size ( 1 ) / 2 ;      % define the center 
Rpic = 0.5 * sqrt ( 2 ) * Size ( 1 ) ; 

tmax = round ( Rpic * 1.1 ) ;
t_int = 0.1 ;
t_range =  -tmax : t_int : tmax ;
t_range = t_range';                                                  % GPU read column vector
Lt = length ( t_range ) ;

%% formualr projection

% R = FormProjectionParal ( Size , thetaRange , t_range , 2 ) ;  % high contrast
% R = reshape ( R' , 1 , Ltheta * Lt ) ;       %  to be consistent with sysmatix

%% compute system matrix

% SysMatrix = GenSysMatParal ( height , width , Size , Center_x , Center_y , thetaRange , t_range ) ;
% picvector = reshape ( flipud( pic )' , height * width , 1  ) ;  % original image
% R = SysMatrix * picvector ;        % generate projection with system matrix

% SysMatrix2 = SysMatrix ;                 % for MART
% SysMatrix2 ( SysMatrix2 ~= 0 ) = 1 ;                 % for MART
% load A.mat
% SysMatrix = A ; 
%% GPU accelerate the projection
picvector = reshape ( flipud( pic ) , height * width , 1  ) ;  % original image
R = ProjectionParallel_2D ( picvector , height , width , Size , thetaRange , t_range ) ; 
R = reshape ( R , Lt , Ltheta ) ;

%% iterative
% S.Kaczmarz Method
% 
% Times = 5 ;
% IterativeTime = 1  ;      % times to iterative
% Lamda = 0.3 ;                       % relaxtion factor
% % Display = zeros ( 1 , height * width ) ;          % store the reconstruction
% Display = ones ( height * width  , 1 ) ;          % store the reconstruction for MART
% MSE = zeros( 1 ,  Times ) ;                                       % judgement parameter
% Picmean = mean2 ( pic ) ;                           % mean of original pic
% 
% while ( IterativeTime <= Times )            % end condition of loop
%     disp ( IterativeTime ) ;
%     for LprojIndex = 1 : Ltheta * Lt           
%             LthetaIndex  = floor ( ( LprojIndex - 1 ) / Lt ) + 1 ; 
%             LtIndex  = mod ( ( LprojIndex - 1 ) , Lt ) + 1 ;
%             Aj = SysMatrix ( LprojIndex , : ) ; 
% %             Bj = SysMatrix2 ( LprojIndex , : ) ;                     % for MART
%             if (  Aj * Aj'  ~= 0 )
%                     Display = Display - Lamda * Aj' .* ( ( Aj * Display ) - R ( ( LthetaIndex - 1 ) * Lt +  LtIndex ) ) ./ ( Aj * Aj' ) ;        
%             end  
% %              if ( Aj * Display' ~= 0 )
% %                     Display = Display +  Display .* Bj * ( R ( LthetaIndex , LtIndex ) - Aj * Display' ) / ( Aj * Display' ) ;     % for MART
% %              end
%              
%     end
% %     Dismean = mean ( Display ) ;
%     MSE ( IterativeTime ) = sum ( sum ( ( Display - picvector ).^2 ) ) /  ( height * width ) ;
%     IterativeTime = IterativeTime + 1 ;
% end
% Display = reshape ( Display , width , height ) ;
% Display =  flipud ( Display' ) ;
% 
% figure , imshow ( Display , [ 0  1 ] ) ;                     % display results
% 
% figure, plot ( 1 : Times , MSE( 1  : Times ) ) ;                          % display error graph
% 
% figure,plot( 1 : size ( pic , 1 ) , Display ( : , 129 ) , 1 : size ( pic , 1 ) , pic ( : , 129 ) ) ;    % display transversal
% title ( ' grey distrubition ' ) ;
% axis ( [ 0 256 0 1 ] ) ;
% 

 toc
 
%%  pixel-drived matrix establish 
% 
 % for thetaindex = 1 : Ltheta
% 
%     thetaRadian =  thetaRange ( thetaindex ) * pi / 180 ;       % angle to radian            
%     
%     for Hindex = 1 : height
% 
%                 for Windex  = 1 : width
%                         Point_x  = ( Windex - 0.5 ) * Resolution - Center_x ;                % pixel central cordinate before rotation 
%                         Point_y  = ( Hindex - 0.5 ) * Resolution - Center_y ; 
%                         t = Point_x * cos ( thetaRadian ) + Point_y * cos ( thetaRadian ) ;      % compute the radial line through the pixel 
%                         t1index = floor ( ( t + tmax ) / t_int ) + 1 ;
%                         t2index = t1index + 1 ;
%                        if ( t >= min ( t_range ) && t < max ( t_range ) ) 
%                                     Ladd1 = ( thetaindex  - 1 ) * Ltheta + t1index ;  Ladd2 = ( thetaindex  - 1 ) * Ltheta + t2index ;      % thetaindex * Ltheta + tindex
%                                     Jadd = ( Windex - 1 ) * width + Hindex ;
%                                     t1 = t_range ( t1index ) ; t2 = t_range ( t2index ) ;  
%                                     if ( thetaRadian >= 0 || thetaRadian < pi / 4 || thetaRadian >= 7 * pi / 4 || thetaRadian < 2 * pi )
%                                         w1 = height / cos ( thetaRadian ) * ( t2 - t ) / t_int ; w2 = height / cos ( thetaRadian ) * ( t - t1 ) / t_int ;
%                                    elseif ( thetaRadian >= 7 * pi / 4 || thetaRadian < 2 * pi )
%                                         w1 = width / cos ( pi / 2 - thetaRadian ) * ( t2 - t ) / t_int ; w2 = width / cos ( pi / 2 - thetaRadian ) * ( t - t1 ) / t_int ;
%                                    elseif ( thetaRadian >= pi / 4 || thetaRadian < 3 * pi / 4 )
%                                         w1 = width / cos ( 3 * pi / 2 - thetaRadian ) * ( t2 - t ) / t_int ; w2 = width / cos ( 3 * pi / 2 - thetaRadian ) * ( t - t1 ) / t_int ;  
%                                    elseif ( thetaRadian >= 5 * pi / 4 || thetaRadian < 7 * pi / 4 )
%                                         w1 = height / cos ( pi - thetaRadian ) * ( t2 - t ) / t_int ; w2 = height / cos (pi - thetaRadian ) * ( t - t1 ) / t_int ;
%                                     end
% %                                    L = [ L ; Ladd1 ; Ladd2 ] ; J = [ J ; Jadd ; Jadd ] ; W = [ W ; w1 ; w2 ] ;
%                                     L ( end + 1 ) = Ladd1 ; L ( end + 1 ) = Ladd2 ; J ( end + 1 ) = Jadd ; J ( end + 1 ) = Jadd ; W ( end + 1 ) = w1 ; W ( end + 1 ) = w2 ; 
%                        end
% 
%                 end%Windex
%     end%Hindex
% 
% end%thetaindex
%% change the pic matrix into row-vector
% pic_row_vector = zeros ( 1 , height * width ) ;
% pic = flipud ( pic ) ;
% for rowindex = 1 : height
%     
%     pic_row_vector ( ( rowindex - 1 ) * width + 1 : rowindex * width ) = pic ( rowindex , : ) ;
%     
% end
%% rebuild the projection matrix to justify the accuracy of system matrix

% Rs = SysMatrix * pic_row_vector' ;
% Rs = reshape ( Rs , Lt , Ltheta ) ;







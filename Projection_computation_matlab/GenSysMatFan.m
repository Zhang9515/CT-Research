% compiled in 2018/01/15 by ZXZ 
% The function generates the system matrix using three vetcors of sparse matrix: index of projection lines , index of pixels , weight
% output : SysMatrix = sparse ( L , J , W , Ltheta * Lt ,  height * width ) ;
% row index of SysMatrix starts with  t1, t2 , t3......tn for theta1 and so on
% column index of SysMatrix starts with x1, x2, x3......xn for y1 and so on 
function [ SysMatrix ] = GenSysMatFan ( height, width, Size, BetaScanRange, Pdomain, RScan, Center_x , Center_y)

    LBeta = length ( BetaScanRange ) ; 
    LP = length ( Pdomain ) ; 
    
    L = [] ; J = [] ; W = [] ;     
    
    Smax = 2 * RScan ; 
    Resolution = max ( Size ) / height ;   % define the resolution of pic
    
    for Betaindex = 1 : LBeta 

                Beta =  BetaScanRange ( Betaindex ) ; 
                source_x = Center_x + RScan * sin( -Beta ) ;
                source_y = Center_y + RScan * cos( Beta ) ;
                
            for Pindex = 1 : LP 
                
                        P = Pdomain ( Pindex ) ;
                        Theta = atan( P / RScan ) ;
                        DetectPoint_xend = Center_x + Smax * sin(Theta) * cos ( -Beta ) + (RScan - Smax * cos(Theta)) * sin(-Beta)  ;      % define end detect point in matlab coordinate
                        DetectPoint_yend = Center_y - Smax * sin(Theta) * sin ( -Beta ) + (RScan - Smax * cos(Theta)) * cos(-Beta)   ;     % Note that : 0 is the start
                        DetectPoint_terminal = [ source_y , source_x ; DetectPoint_yend , DetectPoint_xend ] ;  
                        
                        [ Lx , Ly , x_signal , y_signal , ProjectionLine ] = Siddon2D ( DetectPoint_terminal , Resolution , Size , height , width , RScan ) ;    % call Siddon2D function
                        if ( size( ProjectionLine , 2 ) <=1 )
                            continue ;
                        end
                        X2Y = ( DetectPoint_yend - source_y ) / ( DetectPoint_xend - source_x + eps ) ;
                        Y2X = ( DetectPoint_xend - source_x ) / ( DetectPoint_yend - source_y + eps ) ;
                        X2Y = rangelimit( X2Y ) ; Y2X = rangelimit( Y2X ) ; 
%                         % for the situation of Ly = 0 or Lx = 0 and also
%                         % define the start point in the picture exactly
                        if ( ProjectionLine ( 3 , 1 ) == 1 )
                                DetectPoint_x =  source_x + ( ProjectionLine ( 2 , 1 ) * Resolution - source_y ) * Y2X ;                   % start point of the line
                                DetectPoint_x = ceil ( DetectPoint_x / Resolution ) ;         % to determine the pixel index
                        else
                                DetectPoint_y =  source_y + ( ProjectionLine ( 2 , 1 ) * Resolution - source_x ) * X2Y ;          
                                DetectPoint_y = ceil ( DetectPoint_y / Resolution ) ;        % to determine the pixel index
                        end
%                         
                        for n = 1 : ( Lx + Ly )

                                if ( ProjectionLine ( 3 , n ) == 1 && Ly ~= 0 )                                                              % define the pixel 
                                        DetectPoint_y = ProjectionLine ( 2 , n ) + 0.5 + y_signal * 0.5 ; 
                                elseif ( ProjectionLine ( 3 , n ) == 2 && Lx ~= 0 )
                                        DetectPoint_x = ProjectionLine ( 2 , n ) + 0.5 + x_signal * 0.5 ; 
                                end
%                                 disp(length ( ProjectionLine ))
%                                 disp(n)
%                                 disp('----------------------------------')
%                                 if ( ProjectionLine ( 3 , n ) == 1 )
%                                     D1x =  source_x + ( ProjectionLine ( 2 , n ) * Resolution - source_y ) / tan ( Beta ) ; 
%                                     D1y = ProjectionLine ( 2 , n ) * Resolution ;
%                                 elseif ( ProjectionLine ( 3 , n ) == 2 )
%                                     D1x = ProjectionLine ( 2 , n ) * Resolution ;
%                                     D1y = source_y + ( ProjectionLine ( 2 , n ) * Resolution - source_x ) * tan ( Beta ) ;     
%                                 end
%                                 if ( ProjectionLine ( 3 , n+1 ) == 1 )
%                                     D2x =  source_x + ( ProjectionLine ( 2 , n+1 ) * Resolution - source_y ) / tan ( Beta ) ; 
%                                     D2y = ProjectionLine ( 2 , n+1 ) * Resolution ;
%                                 elseif ( ProjectionLine ( 3 , n+1 ) == 2 )
%                                     D2x = ProjectionLine ( 2 , n+1 ) * Resolution ;
%                                     D2y = source_y + ( ProjectionLine ( 2 , n+1 ) * Resolution - source_x ) * tan ( Beta ) ;     
%                                 end
%                                 
%                                 DetectPoint_x = ceil ( 0.5 * ( D1x + D2x ) / Resolution ) ;
%                                 DetectPoint_y = ceil ( 0.5 * ( D1y + D2y ) / Resolution ) ;
                                
                                if ( DetectPoint_y > 0 && DetectPoint_y <= height && DetectPoint_x > 0 && DetectPoint_x <= width )       
                                         L ( end + 1 ) = ( Betaindex - 1 ) * LP +  Pindex ; J ( end + 1 ) = ( DetectPoint_y - 1 ) * width + DetectPoint_x ;   
                                        W ( end + 1 ) = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                                        % this operation of ( end + 1 ) is more efficient than others 
                                end 
                                
                         end% n 
            end%  j
    end%  thetaindex
    
    SysMatrix = sparse ( L , J , W , LBeta * LP ,  height * width ) ;
    
    
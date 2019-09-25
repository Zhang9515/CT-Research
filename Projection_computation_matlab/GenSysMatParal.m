% compiled in 2018/01/15 by ZXZ 
% The function generates the system matrix using three vetcors of sparse matrix: index of projection lines , index of pixels , weight
% output : SysMatrix = sparse ( L , J , W , Ltheta * Lt ,  height * width ) ;
% row index of SysMatrix starts with  t1, t2 , t3......tn for theta1 and so on
% column index of SysMatrix starts with x1, x2, x3......xn for y1 and so on 
function [ SysMatrix ] = GenSysMatParal ( height , width , Size , Center_x , Center_y ,  thetaRange  , t_range )

    
    
    Ltheta = length ( thetaRange ) ; 
    Lt = length ( t_range ) ; 
    
    L = [] ; J = [] ; W = [] ;     
    
    tmax = t_range ( end ) ;
    Resolution = max ( Size ) / height ;   % define the resolution of pic
    
    for thetaindex = 1 : Ltheta 

                theta =  thetaRange ( thetaindex ) ; 
                
            for tindex = 1 : Lt 
                
                        t = t_range ( tindex ) ;
                        Smax = tmax ; 
                        
                        DetectPoint_xstart = Center_x + t * cos ( -theta ) - Smax * sin ( -theta ) ;    % according to euler equation  
                        DetectPoint_ystart = Center_y - t * sin ( -theta ) - Smax * cos ( -theta ) ;
                        DetectPoint_xend = Center_x + t * cos ( -theta ) + Smax * sin ( -theta ) ;      % define end detect point in matlab coordinate
                        DetectPoint_yend = Center_y - t * sin ( -theta ) + Smax * cos ( -theta ) ;     % Note that : 0 is the start
                        DetectPoint_terminal = [ DetectPoint_ystart , DetectPoint_xstart ; DetectPoint_yend , DetectPoint_xend ] ;  
                        
                        X2Y = ( DetectPoint_yend - DetectPoint_ystart ) / ( DetectPoint_xend - DetectPoint_xstart + eps ) ;
                        Y2X = ( DetectPoint_xend - DetectPoint_xstart ) / ( DetectPoint_yend - DetectPoint_ystart + eps ) ;
                        X2Y = rangelimit( X2Y ) ; Y2X = rangelimit( Y2X ) ; 
                        
                        [ Lx , Ly , x_signal , y_signal , ProjectionLine ] = Siddon2D ( DetectPoint_terminal , Resolution , Size , height , width , Smax ) ;    % call Siddon2D function
                        
                        % for the situation of Ly = 0 or Lx = 0 and also
                        % define the start point in the picture exactly
                        if ( ProjectionLine ( 3 , 1 ) == 1 )
                                DetectPoint_x =  DetectPoint_xstart + ( ProjectionLine ( 2 , 1 ) * Resolution - DetectPoint_ystart ) * Y2X ;                   % start point of the line
                                DetectPoint_x = ceil ( DetectPoint_x / Resolution ) ;         % to determine the pixel index
                        else
                                DetectPoint_y =  DetectPoint_ystart + ( ProjectionLine ( 2 , 1 ) * Resolution - DetectPoint_xstart ) * X2Y ;          
                                DetectPoint_y = ceil ( DetectPoint_y / Resolution ) ;        % to determine the pixel index
                        end
                        
                        for n = 1 : length ( ProjectionLine )

                                if ( ProjectionLine ( 3 , n ) == 1 && Ly ~= 0 )                                                              % define the pixel 
                                        DetectPoint_y = ProjectionLine ( 2 , n ) + 0.5 + y_signal * 0.5 ; 
                                elseif ( ProjectionLine ( 3 , n ) == 2 && Lx ~= 0 )
                                        DetectPoint_x = ProjectionLine ( 2 , n ) + 0.5 + x_signal * 0.5 ; 
                                end

                                if ( DetectPoint_y > 0 && DetectPoint_y <= height && DetectPoint_x > 0 && DetectPoint_x <= width )       
                                         L ( end + 1 ) = ( thetaindex - 1 ) * Lt +  tindex ; J ( end + 1 ) = ( DetectPoint_y - 1 ) * width + DetectPoint_x ;   
                                        W ( end + 1 ) = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                                        % this operation of ( end + 1 ) is more efficient than others 
                                end 
                                
                         end% n 
            end%  j
    end%  thetaindex
    
    SysMatrix = sparse ( L , J , W , Ltheta * Lt ,  height * width ) ;
    
    
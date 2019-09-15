% compiled in 2018/01/15 by ZXZ
% grid-idea  Siddon
function [ Lx , Ly , x_signal , y_signal , ProjectionLine ] = Siddon2D ( DetectPoint_terminal , Resolution , Size , height , width , Smax )
% Resolution denotes length of each unit 
% start point denotes tvertical of the smallest value,  Note that : 0 is the start
% L , J , W are three vetcors of sparse matrix: index of projection lines , index of pixels , weight, repectively
% 
% 
%      tv \      |y  /t
%             \   | /
% -----------|/----------->x
%                  |

DetectPoint_xstart = DetectPoint_terminal ( 1 , 2 ) ; DetectPoint_ystart = DetectPoint_terminal ( 1 , 1 ) ; 
DetectPoint_xend = DetectPoint_terminal ( 2 , 2 ) ; DetectPoint_yend = DetectPoint_terminal ( 2 , 1 ) ; 

xstartindexL = floor ( DetectPoint_xstart / Resolution ) ; xstartindexR = ceil ( DetectPoint_xstart / Resolution ) ;
ystartindexL = floor ( DetectPoint_ystart / Resolution ) ; ystartindexR = ceil ( DetectPoint_ystart / Resolution ) ;
xendindexL = floor ( DetectPoint_xend / Resolution ) ; xendindexR = ceil ( DetectPoint_xend / Resolution ) ;
yendindexL = floor ( DetectPoint_yend / Resolution ) ; yendindexR = ceil ( DetectPoint_yend / Resolution ) ;

% to determine the range of y
    if ( yendindexR == ystartindexR || yendindexL == ystartindexL )    
            y_range = [ ] ;
            y_signal = 0 ;
    else
            if ( DetectPoint_yend > DetectPoint_ystart )
                     y_signal = 1 ;   
                     if ( DetectPoint_ystart < 0 )
                             ystart = 0 ; 
                     else
                             ystart = ystartindexL ;    % if ystart is within  y range of the picture, that means the projection line will first cross the x-grid  
                     end
                     if ( DetectPoint_yend > Size ( 1 ) )
                             yend = height ;
                     else
                             yend = yendindexR ;
                     end                         
            elseif ( DetectPoint_yend < DetectPoint_ystart )
                    y_signal = -1 ; 
                    if ( DetectPoint_yend < 0 )
                             ystart = 0 ;
                     else
                             ystart = yendindexL ;
                     end
                     if ( DetectPoint_ystart > Size ( 1 ) )
                             yend = height ;
                     else
                             yend = ystartindexR ;
                     end
            end
            y_range = ystart : yend ;     % we should always keep end > start
    end
    % to determine the range of x
    if ( xendindexR == xstartindexR || xendindexL == xstartindexL )
             x_range = [ ] ;
             x_signal = 0 ;
    else
            if ( DetectPoint_xend > DetectPoint_xstart )
                    x_signal = 1 ; 
                    if ( DetectPoint_xstart < 0 )
                             xstart = 0 ;
                     else
                             xstart = floor ( DetectPoint_xstart / Resolution ) ;
                     end
                     if ( DetectPoint_xend > Size ( 2 ) )
                             xend = width ;
                     else
                             xend = ceil ( DetectPoint_xend / Resolution ) ;
                     end
            elseif ( DetectPoint_xend < DetectPoint_xstart )
                    x_signal = -1 ; 
                    if ( DetectPoint_xend < 0 )
                             xstart = 0 ;
                     else
                             xstart = floor ( DetectPoint_xend / Resolution ) ;
                     end
                     if ( DetectPoint_xstart > Size ( 2 ) )
                             xend = width ;
                     else
                             xend = ceil ( DetectPoint_xstart / Resolution ) ;
                     end
            end
            x_range = xstart : xend ;              % we should always keep end > start
    end

    Ly = length ( y_range ) ; Lx = length ( x_range ) ;   
    ProjectionLine = zeros ( 3 , Ly  + Lx ) ;
    % record the proportion of projection line decided by cross point
    ProjectionLine ( 1 , 1 : Ly  ) = abs ( ( y_range * Resolution - DetectPoint_ystart ) / ( DetectPoint_yend - DetectPoint_ystart ) ) * 2 * Smax ;
    ProjectionLine ( 2 , 1 : Ly  ) = y_range ;
    ProjectionLine ( 3 , 1 : Ly  ) = 1 ;                              % 1 represents y
    ProjectionLine ( 1 , ( Ly  + 1 ) : ( Ly  + Lx ) ) = abs ( ( x_range * Resolution - DetectPoint_xstart ) / ( DetectPoint_xend - DetectPoint_xstart ) ) * 2 * Smax ;
    ProjectionLine ( 2 , ( Ly  + 1 ) : ( Ly  + Lx ) ) = x_range ;
    ProjectionLine ( 3 , ( Ly  + 1 ) : ( Ly  + Lx ) ) = 2 ;                              % 2 represents x                   
    ProjectionLine = sortrows ( ProjectionLine' )' ;

    
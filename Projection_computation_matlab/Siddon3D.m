% compiled in 18/03/12 by ZXZ 
% grid-idea  Siddon algrithim
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

function [ Lt , Ls , Lz , t_signal , s_signal , z_signal , ProjectionLine ] = Siddon3D ( DetectPoint_terminal , Resolution , Size , Smax ) 

                        source_t = DetectPoint_terminal ( 1 , 1 ) ; source_s = DetectPoint_terminal ( 1 , 2 ) ; source_z = DetectPoint_terminal ( 1 , 3 ) ;  
                        DetectPoint_tend = DetectPoint_terminal ( 2 , 1 ) ; DetectPoint_send = DetectPoint_terminal ( 2 , 2 ) ; DetectPoint_zend = DetectPoint_terminal ( 2 , 3 ) ;
                        
                        % record the t,s,z grids on two sides of start/end points 
                        source_tindexL = floor ( source_t / Resolution(1) ) ; source_tindexR = ceil ( source_t / Resolution(1) ) ;
                        source_sindexL = floor ( source_s / Resolution(2) ) ; source_sindexR = ceil ( source_s / Resolution(2) ) ;
                        source_zindexL = floor ( source_z / Resolution(3) ) ; source_zindexR = ceil ( source_z / Resolution(3) ) ;
                        tendindexL = floor ( DetectPoint_tend / Resolution(1) ) ; tendindexR = ceil ( DetectPoint_tend / Resolution(1) ) ;
                        sendindexL = floor ( DetectPoint_send / Resolution(2) ) ; sendindexR = ceil ( DetectPoint_send / Resolution(2) ) ;
                        zendindexL = floor ( DetectPoint_zend / Resolution(3) ) ; zendindexR = ceil ( DetectPoint_zend / Resolution(3) ) ;

                        % to determine the range of t
                        if ( tendindexR == source_tindexR || tendindexL == source_tindexL )    
                                t_range = [ ] ;
                                t_signal = 0 ;
                        else
                                if ( DetectPoint_tend > source_t )
                                         t_signal = 1 ;   
                                         if ( source_t < 0 )
                                                 tstart = 0 ; 
                                         else
                                                 tstart = source_tindexL ;    % if tstart is within  t range of the picture, that means the projection line will first cross the x-grid  
                                         end
                                         if ( DetectPoint_tend > Size (1) * Resolution(1) )
                                                 tend = Size (1) ;
                                         else
                                                 tend = tendindexR ;
                                         end                         
                                elseif ( DetectPoint_tend < source_t )
                                        t_signal = -1 ; 
                                        if ( DetectPoint_tend < 0 )
                                                 tstart = 0 ;
                                         else
                                                 tstart = tendindexL ;
                                         end
                                         if ( source_t > Size (1)  * Resolution(1) )
                                                 tend = Size (1) ;
                                         else
                                                 tend = source_tindexR ;
                                         end
                                end
                                t_range = tstart : tend ;     % we should always keep end > start
                        end
                        % to determine the range of s
                        if ( sendindexR == source_sindexR || sendindexL == source_sindexL )    
                                s_range = [ ] ;
                                s_signal = 0 ;
                        else
                                if ( DetectPoint_send > source_s )
                                         s_signal = 1 ;   
                                         if ( source_s < 0 )
                                                 sstart = 0 ; 
                                         else
                                                 sstart = source_sindexL ;    % if sstart is within  s range of the picture, that means the projection line will first cross the x-grid  
                                         end
                                         if ( DetectPoint_tend > Size (2)  * Resolution(2) )
                                                 send = Size (2) ;
                                         else
                                                 send = sendindexR ;
                                         end                         
                                elseif ( DetectPoint_send < source_s )
                                        s_signal = -1 ; 
                                        if ( DetectPoint_send < 0 )
                                                 sstart = 0 ;
                                         else
                                                 sstart = sendindexL ;
                                         end
                                         if ( source_s > Size (2)  * Resolution(2) )
                                                 send = Size (2) ;
                                         else
                                                 send = source_sindexR ;
                                         end
                                end
                                s_range = sstart : send ;     % we should always keep end > start
                        end
                        % to determine the range of z
                        if ( zendindexR == source_zindexR || zendindexL == source_zindexL )    
                                z_range = [ ] ;
                                z_signal = 0 ;
                        else
                                if ( DetectPoint_zend > source_z )
                                         z_signal = 1 ;   
                                         if ( source_z < 0 )
                                                 zstart = 0 ; 
                                         else
                                                 zstart = source_zindexL ;    % if ystart is within  y range of the picture, that means the projection line will first cross the x-grid  
                                         end
                                         if ( DetectPoint_zend > Size (3)  * Resolution(3) )
                                                 zend = Size (3) ;
                                         else
                                                 zend = zendindexR ;
                                         end                         
                                elseif ( DetectPoint_zend < source_z )
                                        z_signal = -1 ; 
                                        if ( DetectPoint_zend < 0 )
                                                 zstart = 0 ;
                                         else
                                                 zstart = zendindexL ;
                                         end
                                         if ( source_z > Size (3)  * Resolution(3) )
                                                 zend = Size (3) ;
                                         else
                                                 zend = source_zindexR ;
                                         end
                                end
                                z_range = zstart : zend ;     % we should always keep end > start
                        end
                        
                        Lt = length ( t_range ) ; Ls = length ( s_range ) ; Lz = length ( z_range ) ; 
                        ProjectionLine = zeros ( 3 , Lt + Ls + Lz ) ;
                        ProjectionLine ( 1 , 1 : Lt ) = abs ( ( t_range * Resolution(1) - source_t ) / ( DetectPoint_tend - source_t + 1e-10 ) ) * Smax ;
                        ProjectionLine ( 2 , 1 : Lt ) = t_range ;
                        ProjectionLine ( 3 , 1 : Lt ) = 1 ;                              % 1 represents t
                        ProjectionLine ( 1 , ( Lt + 1 ) : ( Ls + Lt ) ) = abs ( ( s_range * Resolution(2) - source_s ) / ( DetectPoint_send - source_s + 1e-10 ) ) * Smax ;
                        ProjectionLine ( 2 , ( Lt + 1 ) : ( Ls + Lt ) ) = s_range ;
                        ProjectionLine ( 3 , ( Lt + 1 ) : ( Ls + Lt ) ) = 2 ;                              % 2 represents s
                        ProjectionLine ( 1 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = abs ( ( z_range * Resolution(3) - source_z ) / ( DetectPoint_zend - source_z + 1e-10 ) ) * Smax ;
                        ProjectionLine ( 2 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = z_range ;
                        ProjectionLine ( 3 , ( Lt + Ls + 1 ) : ( Ls + Lt + Lz ) ) = 3 ;                              % 3 represents z
                        ProjectionLine = sortrows ( ProjectionLine' )' ;
                        

clear all;
% close all;
%%
% parameter define
theta_0 = 1 ;                                     %角度间隔
N = 360 / theta_0 ;                                    %角度探测次数
% pic = phantom ( 512 ) ;
pic = ones ( 512 , 512 ) ;                                  %内置图形测试
% figure , imshow ( pic ) ;                                %以伪彩色形式显示
%figure,imshow(p,[]); 
title ( ' Original image ' ) ;
theta = theta_0 : theta_0 : 360 ;                                %radon scanning range
[ height , width ] = size ( pic ) ;              % store the size of picture
Center_i = height / 2 ;  Center_j = width / 2 ;      % define the center 
tmax = round ( 0.5 * sqrt ( 2 ) * size ( pic , 1 ) * 1.1 ) ;
% tmax = 364 ; 
t_int = 1 ;
t_range =  -tmax : t_int : tmax ;

%% my radon

R = zeros ( length ( t_range ) ,  length ( theta ) ) ;   % create space to store fan projection
Smax = 2 * tmax ;
for theta_axis = 1 : length ( theta )
%         theta_axis = 90 ;
        
        thetaradian = theta ( theta_axis ) * pi / 180 ;
        Htransform ( 1 , 1 ) = cos ( thetaradian ) ;
        Htransform ( 1 , 2 ) = -sin ( thetaradian ) ;
        Htransform ( 2 , 1 ) = sin ( thetaradian ) ;
        Htransform ( 2 , 2 ) = cos ( thetaradian ) ;
        for t_axis = 1 : length ( t_range )
%                 t_axis = 109 ; 
                T = t_range ( t_axis ) ;
                Pstart = [ -tmax ; T ] ;
                Pend = [ tmax ; T ] ; 
                Resultstart = Htransform * Pstart ;
                Resultend = Htransform * Pend ;
                dect_istart = Resultstart ( 1 ) + Center_i ;
                dect_jstart = Resultstart ( 2 ) + Center_j ;
                dect_iend = Resultend ( 1 ) + Center_i ;
                dect_jend = Resultend ( 2 ) + Center_j ;
                
                % define projection line direction and range of i and j
                if ( floor ( dect_iend ) == floor ( dect_istart ) || ceil ( dect_iend ) == ceil ( dect_istart ) )    
                    i_range = [ ] ;
                elseif ( dect_iend > dect_istart )                                                         
                     i_signal = 1 ;   
                     i_range = ceil ( dect_istart ) : floor ( dect_iend ) ;
                elseif ( dect_iend < dect_istart )
                    i_signal = -1 ; 
                    i_range = ceil ( dect_iend ) : floor ( dect_istart ) ;
                end
                if ( floor ( dect_jend ) == floor ( dect_jstart ) || ceil ( dect_jend ) == ceil ( dect_jstart ) )
                     j_range = [ ] ;
                elseif ( dect_jend > dect_jstart )
                    j_signal = 1 ; 
                    j_range = ceil ( dect_jstart ) : floor ( dect_jend ) ;
                elseif ( dect_jend < dect_jstart )
                    j_signal = -1 ; 
                    j_range = ceil ( dect_jend ) : floor ( dect_jstart ) ;
                end       
                
                Li = length ( i_range ) ; Lj = length ( j_range ) ;
                ProjectionLine = zeros ( 3 , Li + Lj ) ;
                ProjectionLine ( 1 , 1 : Li ) = abs ( ( i_range - dect_istart ) / ( dect_iend - dect_istart ) ) * Smax ;
                ProjectionLine ( 2 , 1 : Li ) = i_range ;
                ProjectionLine ( 3 , 1 : Li ) = 1 ;                              % 1 represents i
                ProjectionLine ( 1 , ( Li + 1 ) : ( Li + Lj ) ) = abs ( ( j_range - dect_jstart ) / ( dect_jend - dect_jstart ) ) * Smax ;
                ProjectionLine ( 2 , ( Li + 1 ) : ( Li + Lj ) ) = j_range ;
                ProjectionLine ( 3 , ( Li + 1 ) : ( Li + Lj ) ) = 2 ;                              % 2 represents j
                ProjectionLine = sortrows ( ProjectionLine' )' ;
                DetectPoint_i = ceil ( dect_istart ) ; DetectPoint_j = ceil ( dect_jstart ) ;                   % start point of the line        

                for n = 1 : length ( ProjectionLine )                                                                                    
                    if ( ProjectionLine ( 3 , n ) == 1 && Li ~= 0 )                                                               % define the pixel 
                            DetectPoint_i = ProjectionLine ( 2 , n ) + 0.5 + i_signal * 0.5 ; 
                    elseif ( ProjectionLine ( 3 , n ) == 2 && Lj ~= 0 )
                            DetectPoint_j = ProjectionLine ( 2 , n ) + 0.5 + j_signal * 0.5 ; 
                    end
                    
                    if ( Li == 0 && dect_istart == 0 && DetectPoint_j > 0 && DetectPoint_j <= width )                % boundary condition
                        f = pic ( 1 , DetectPoint_j ) / 2 ;
                        LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end   
                    if ( Lj == 0 && dect_jstart == 0 && DetectPoint_i > 0 && DetectPoint_i <= width )
                        f = pic ( DetectPoint_i , 1 ) / 2 ;
                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                        R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end   
                    
                    if ( DetectPoint_i > 0 && DetectPoint_i <= height && DetectPoint_j > 0 && DetectPoint_j <= width ) 
                         
                        if ( Li == 0 && dect_istart == DetectPoint_i )                        % condition when projection line overlay the axis line                
                             if ( dect_istart == height )
                                 f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i + 1 , DetectPoint_j ) ) / 2 ;
                             end
                         elseif ( Lj == 0 && dect_jstart == DetectPoint_j ) 
                             if ( dect_jstart == width )
                                 f = pic ( DetectPoint_i , DetectPoint_j ) / 2 ;
                             else 
                                 f = ( pic ( DetectPoint_i , DetectPoint_j ) + pic ( DetectPoint_i , DetectPoint_j + 1 ) ) / 2 ;
                             end
                         else 
                                f = pic ( DetectPoint_i , DetectPoint_j ) ;
                         end
                       
                         LengthUnit = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;                     
                         R ( t_axis , theta_axis ) = R ( t_axis , theta_axis ) + f * LengthUnit ;    
                    end
                end
        end
end
figure , imshow ( R ) ;    
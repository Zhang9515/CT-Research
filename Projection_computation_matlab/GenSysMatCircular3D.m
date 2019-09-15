% compiled by ZXZ in 2018/03/13
% flat-panel
% Size = [ Lt , Ls , Lz ]
% row index of SysMatrix starts with  P1, P2 , P3......Pn for Xigama1 and
%      Xigama2 ... Xigamam for Beta1 and so on
% column index  of SysMatrix starts with t1, t2, t3......tn for s1 and s2
%    ... sm for z1 and so on
function [ SysMatrix ] = GenSysMatCircular3D ( Size , Center_t , Center_s , Center_z , Pdomain  , Xigamadomain , BetaScanRange,  Distance , Resolution )

LBeta = length ( BetaScanRange ) ;
LXigama = length ( Xigamadomain ) ;
LP = length ( Pdomain ) ;

L = [] ; J = [] ; W = [] ;     

for Betaindex = 1 : LBeta 
            for Xigamaindex = 1 : LXigama
                    for Pindex = 1 : LP
                        
                            beta =  BetaScanRange ( Betaindex ) ;
                            betaRadian = beta * pi / 180 ;                                          % angle to radian
                            source_t = Center_t - Distance * sin ( betaRadian ) ;      %define the source in matlab coordinate
                            source_s = Center_s + Distance * cos ( betaRadian ) ;  
                            source_z  = Center_z ;                                                        

                                gama = atan ( Xigamadomain ( Xigamaindex ) / Distance ) ;          % radian angle in s-z coordinate plane
                                Distance_shift = Distance / cos ( gama ) ;                    % length of DO'

                                    theta = atan ( Pdomain ( Pindex ) / Distance_shift ) ;                                    % radian angle in s'-t coordinate plane 

                                    Smax = 2 * Distance ;         % define the maximum length of the projection line
                                    DetectPoint_tend = Center_t + Smax * sin ( theta ) * cos ( betaRadian ) - ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * sin ( betaRadian ) ; 
                                    DetectPoint_send = Center_s + Smax * sin ( theta ) * sin ( betaRadian ) + ( Distance - Smax * cos ( theta ) * cos ( gama ) ) * cos ( betaRadian ) ;
                                    DetectPoint_zend = Center_z + Smax * cos ( theta ) * sin ( gama ) ;

                                    % record both terminals of the current projection line
                                    DetectPoint_terminal = [ source_t , source_s , source_z ; DetectPoint_tend , DetectPoint_send , DetectPoint_zend ] ;

                                    % call Siddon3D function
                                      [ Lt , Ls , Lz , t_signal , s_signal , z_signal , ProjectionLine ] = Siddon3D ( DetectPoint_terminal , Resolution , Size , Smax ) ;    

                                       % for the situation of Ly = 0 or Lx = 0 or Lz =0 and also
                                        % define the start point in the picture exactly

                                        if ( ProjectionLine ( 3 , 1 ) == 1 )
                                                 % distance between this point and source can be computed since t is known
                                                S_t = ( ProjectionLine ( 2 , 1 ) - Center_t + Distance * sin ( betaRadian ) ) / ( sin ( theta ) * cos ( betaRadian ) + cos ( theta ) * cos ( gama ) * sin ( betaRadian ) ) ;

                                                DetectPoint_s =  Center_s + S_t * sin ( theta ) * sin ( betaRadian ) + ( Distance - S_t * cos ( theta ) * cos ( gama ) ) * cos ( betaRadian ) ;                   % start point of the line
                                                DetectPoint_s = ceil ( DetectPoint_s / Resolution(2) ) ;         % to determine the pixel index
                                                DetectPoint_z =  Center_z + S_t * cos ( theta ) * sin ( gama ) ;             % start point of the line
                                                DetectPoint_z = ceil ( DetectPoint_z / Resolution(3) ) ;         % to determine the pixel index

                                        elseif ( ProjectionLine ( 3 , 1 ) == 2 )
                                                % distance between this point and source can be computed since s is known
                                                S_s = ( ProjectionLine ( 2 , 1 ) - Center_s - Distance * cos ( betaRadian ) ) / ( sin ( theta ) * cos ( betaRadian ) - cos ( theta ) * cos ( gama ) * sin ( betaRadian ) ) ;

                                                DetectPoint_t = Center_t + S_s * sin ( theta ) * cos ( betaRadian ) - ( Distance - S_s * cos ( theta ) * cos ( gama ) ) * sin ( betaRadian ) ; 
                                                DetectPoint_t = ceil ( DetectPoint_t / Resolution(1) ) ;
                                                DetectPoint_z =  Center_z + S_s * cos ( theta ) * sin ( gama ) ;             % start point of the line
                                                DetectPoint_z = ceil ( DetectPoint_z / Resolution(3) ) ;         % to determine the pixel index
                                        elseif ( ProjectionLine ( 3 , 1 ) == 3 ) 
                                                % distance between this point and source can be computed since z is known
                                                S_z = ( ProjectionLine ( 2 , 1 ) - Center_z ) / ( cos ( theta ) * sin ( gama ) ) ;

                                                DetectPoint_s = Center_s + S_z * sin ( theta ) * sin ( betaRadian ) + ( Distance - S_z * cos ( theta ) * cos ( gama ) ) * cos ( betaRadian ) ;                   % start point of the line
                                                DetectPoint_s = ceil ( DetectPoint_s / Resolution(2) ) ;         % to determine the pixel index
                                                DetectPoint_t = Center_t + S_z * sin ( theta ) * cos ( betaRadian ) - ( Distance - S_z * cos ( theta ) * cos ( gama ) ) * sin ( betaRadian ) ; 
                                                DetectPoint_t = ceil ( DetectPoint_t / Resolution(1) ) ;
                                        end

                                        for n = 1 : length ( ProjectionLine )

                                                if ( ProjectionLine ( 3 , n ) == 1 && Lt ~= 0 )                                                              % define the pixel 
                                                        DetectPoint_t = ProjectionLine ( 2 , n ) + 0.5 + t_signal * 0.5 ; 
                                                elseif ( ProjectionLine ( 3 , n ) == 2 && Ls ~= 0 )
                                                        DetectPoint_s = ProjectionLine ( 2 , n ) + 0.5 + s_signal * 0.5 ; 
                                                elseif ( ProjectionLine ( 3 , n ) == 3 && Lz ~= 0 )
                                                        DetectPoint_z = ProjectionLine ( 2 , n ) + 0.5 + z_signal * 0.5 ;
                                                end
                                                
                                                if ( DetectPoint_t > 0 && DetectPoint_t <= Size (1) && DetectPoint_s > 0 && DetectPoint_s <= Size (2) && DetectPoint_z > 0 && DetectPoint_z <= Size (3) )       
                                                     if(DetectPoint_z==59 &&  DetectPoint_t  == 58 && DetectPoint_s == 64 )
                                                        a=0;
                                                     end     
                                                     L ( end + 1 ) = ( Betaindex - 1 ) * LXigama * LP + ( Xigamaindex - 1 ) * LP + Pindex ; J ( end + 1 ) = ( DetectPoint_z - 1 ) * Size(1) * Size(2) + ( DetectPoint_s - 1 ) * Size(1) + DetectPoint_t ;   
                                                        if (ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) > 0 )
                                                            W ( end + 1 ) = ProjectionLine ( 1 , n + 1 ) - ProjectionLine ( 1 , n ) ;  
                                                        else
                                                            W ( end + 1 ) = ProjectionLine ( 1 , n + 2 ) - ProjectionLine ( 1 , n ) ;
                                                        end
                                                        % this operation of ( end + 1 ) is more efficient than others 
                                                end 
                                
                                        end% n 
                    end%P
            end%Xigama
end%Beta

SysMatrix = sparse ( L , J , W , LBeta * LXigama * LP ,  prod ( Size ) ) ;

function H_omega = HstarOmega( omega , ps )
% 2018/11/12
% assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
% X : from left to right; Y : from down to up; Z : from near to far
%input: omega: 3-3-Lx-Ly-Lz, ps: length of patch; output: H_omega: N(Lx*Ly*Lz)-1 
    [~, ~, Lx, Ly, Lz] = size ( omega ) ;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    H = HessianLOGtemplate3D( ps ) ;
%     omega_stack = repmat ( omega , 1, 1, 1, ps, ps, ps) ;
    % 
    H_omega = zeros( Lx, Ly, Lz ) ; 

    % actually we are computing the adjoint of hessian3D matrix
    H_transpose = H ;
    
    % exchange the corresponding row/column between the specific axis
        num = ceil (( ps - 1 ) / 2) ;
        c = H_transpose( : , : , 1 : num , : , : ) ; 
        H_transpose( : , : , 1 : num , : , : ) = H_transpose( : , : , ps - num +1 , : , : ) ; 
        H_transpose( : , : , ps - num +1 , : , : ) = c ;
        
        c = H_transpose( : , : , : , 1 : num , : ) ;
        H_transpose( : , : , : , 1 : num , : ) = H_transpose( : , : , : , ps - num +1 , : ) ; 
        H_transpose( : , : , : , ps - num +1 , : ) = c ; 
        
        c = H_transpose( : , : , : , : , 1 : num ) ;
        H_transpose( : , : , : , : , 1 : num ) = H_transpose( : , : , : , : , ps - num +1 ) ; 
        H_transpose( : , : , : , : , ps - num +1 ) = c ;

%     for UJindex = 1 : sizeofData
%                 H_omega( UJindex ) = sum(sum(sum(sum(sum ( squeeze( omega_stack( : , : , UJindex , : , : , : ) ) .* H_transpose ) ) ) ) ) ; 
%     end       
        for i = 1 : ps
            for j = 1 : ps
                    H_omega = H_omega + convn( omega ( i , j , : , : , : ) , H_transpose( i , j , : , : , : ) , 'same' )  ;               
            end        
        end
        
        H_omega = reshape ( H_omega , Lx * Ly * Lz , 1) ; 


end

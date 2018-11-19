function H_omega = HstarOmega( omega , ps )
% 2018/11/12
% assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
% X : from left to right; Y : from down to up; Z : from near to far
%input: omega: 9*N, ps: length of patch; output: H_omega: N*1 
    [~, sizeofData] = size ( omega ) ;
    H = HessianLOGtemplate3D( ps ) ;
    
    omega_stack = repmat ( omega , 1, 1, ps, ps, ps) ;
    H_omega = zeros( sizeofData , 1 ) ; 

    % actually we are computing the adjoint of hessian3D matrix
    H_transpose = H ;
    
    % exchange the corresponding row/column between the specific axis
    for num = 1 : ( ps - 1 ) / 2
        c = H_transpose( : , [num,3] , : , : ) ; 
        H_transpose( : , [num,3] , : , : ) = H_transpose( : , [num,1] , : , : ) ; 
        H_transpose( : , [num,1] , : , : ) = c ;
        c = H_transpose( : , : , [num,3] , : ) ;
        H_transpose( : , : , [num,3] , : ) = H_transpose( : , : , [num,1] , : ) ; 
        H_transpose( : , : , [num,1] , : ) = c ; 
        c = H_transpose( : , : , : , [num,3] ) ;
        H_transpose( : , : , : , [num,3] ) = H_transpose( : , : , : , [num,1] ) ; 
        H_transpose( : , : , : , [num,1] ) = c ;
    end

    for UJindex = 1 : sizeofData
                H_omega( UJindex ) = sum(sum(sum(sum ( squeeze( omega_stack( : , UJindex , : , : , : ) ) .* H_transpose ) ) ) ) ; 
    end       

end

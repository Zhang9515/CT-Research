function img = cgls4TV ( SysMatrix, divergence, b_CG , iter_CG , miu , lamda, img_previous)

    % because the parameter matrix is miu*A'*A + lamda*(gx'*gx+gy'*gy), so
    % it is symetric and definite. To reduce the computation, here omit the
    % transpose operation on the parameter matrix. Since A is super large,
    % here i split the matrix, letting each component multiply separately.
    times = 1 ; 

    residual = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous ) ; 
    r0 = matrixMultiply( SysMatrix , divergence , miu , lamda , residual ) ;   % initial residual for conjugated gradient algorithm ( CG )
    d0 = r0 ;              % initial search direction for conjugated gradient algorithm ( CG )

    while ( times <= iter_CG )

        
        % conjugated gradient algorithm ( CG )

        if  ( times == 1 )
                d = d0 ;  
                r_previous = r0 ; 
        else
                r_previous = r_next ; 
        end
        AG = matrixMultiply( SysMatrix , divergence , miu , lamda , d ) ; 
        alpha = ( d' * r_previous ) / ( AG' * AG ) ;
        img = img_previous + alpha * d ;
        if ( mod ( times , 10 ) == 1 )
            residual = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img ) ; 
            r_next = matrixMultiply( SysMatrix , divergence , miu , lamda , residual ) ;
        else 
            r_next = r_previous - alpha * matrixMultiply( SysMatrix , divergence , miu , lamda , AG ) ;      
        end 
        beta = ( r_next' * r_next ) / ( r_previous' * r_previous ) ;  
        d = r_next + beta * d ; 
        
        local_e = norm(img-img_previous,2);
        disp ( ['CG iteration: ', num2str(times),'; CG_local_e: ', num2str(local_e)] ) ;
        times = times + 1 ; 
        img_previous = img ; 
    end
end

function output = matrixMultiply( SysMatrix , divergence , miu , lamda , colvect )
    output = miu * (SysMatrix') * ( SysMatrix * colvect ) + lamda * divergence * colvect ;
end



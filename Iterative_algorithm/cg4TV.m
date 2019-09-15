function img = cg4TV ( SysMatrix, divergence, b_CG , iter_CG , miu , lamda, img_previous_iterative)
% using previous result to compute the search direction      
%     
%     img_previous = zeros(size(SysMatrix,2),1) ;     % initialization with previous result
%     b_CG = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous_iterative ) ; 
    img_previous = img_previous_iterative ;
    
    % because the parameter matrix is miu*A'*A + lamda*(gx'*gx+gy'*gy), so
    % it is symetric and definite. To reduce the computation, here omit the
    % transpose operation on the parameter matrix. Since A is super large,
    % here i split the matrix, letting each component multiply separately.

    times = 1 ; 
    threshold = 1e-5 ;
    r0 = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous ) ;   % initial residual for conjugated gradient algorithm ( CG ), which is negative of the gradient
    d0 = r0 ;              % initial search direction for conjugated gradient algorithm ( CG )
    local_e = 100 ;   % initial

    while ( times <= iter_CG && local_e >= threshold )
  
        % conjugated gradient algorithm ( CG )

        if  ( times == 1 )
                d = d0 ;  
                r_previous = r0 ; 
        else
                r_previous = r_next ; 
        end
        Ad = matrixMultiply( SysMatrix , divergence , miu , lamda , d ) ; 
        alpha = ( d' * r_previous ) / ( d' * Ad ) ;
        img = img_previous + alpha * d ;
        residual = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img ) ; 
        if ( mod ( times , 5 ) == 1 )
            r_next = residual ;
        else 
            r_next = r_previous - alpha * Ad ;      
        end 
        beta = ( r_next' * r_next ) / ( r_previous' * r_previous ) ;  
        d = r_next + beta * d ; 
        
        local_e = LocalError(img,img_previous);        
        loss = norm( residual , 2) / 2;   % objective function
        disp ( ['CG iteration: ', num2str(times),';   |    CG_local_e: ', num2str(local_e), ';   |    CG_Loss: ', num2str(loss)] ) ;
        times = times + 1 ; 
        img_previous = img ; 
    end
%     img = img + img_previous_iterative ;
end

function output = matrixMultiply( SysMatrix , divergence , miu , lamda , colvect )
    output = miu * (SysMatrix') * ( SysMatrix * colvect ) + lamda * divergence * colvect ;
end



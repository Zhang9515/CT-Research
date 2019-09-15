% function img = cgls4TV ( SysMatrix, divergence, b_CG , iter_CG , miu , lamda, img_previous_iterative)
% % using previous result to compute the search direction      
% %     
% %     img_previous = zeros(size(SysMatrix,2),1) ;     % initialization with previous result
% %     b_CG = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous_iterative ) ; 
%     img_previous = img_previous_iterative ;
%     
%     % because the parameter matrix is miu*A'*A + lamda*(gx'*gx+gy'*gy), so
%     % it is symetric and definite. To reduce the computation, here omit the
%     % transpose operation on the parameter matrix. Since A is super large,
%     % here i split the matrix, letting each component multiply separately.
% 
%     times = 1 ; 
%     threshold = 1e-5 ;
%     residual = b_CG ; 
%     r0 = matrixMultiply( SysMatrix , divergence , miu , lamda , residual ) ;   % initial gradient for conjugated gradient algorithm ( CG )
%     d0 = r0 ;              % initial search direction for conjugated gradient algorithm ( CG )
%     local_e = 100 ;   % initial
% 
%     while ( times <= iter_CG && local_e >= threshold )
%   
%         % conjugated gradient algorithm ( CG )
% 
%         if  ( times == 1 )
%                 d = d0 ;  
%                 r_previous = r0 ; 
%         else
%                 r_previous = r_next ; 
%         end
%         AG = matrixMultiply( SysMatrix , divergence , miu , lamda , d ) ; 
%         alpha = ( d' * r_previous ) / ( AG' * AG ) ;
%         img = img_previous + alpha * d ;
%         residual = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img ) ; 
%         r_next = matrixMultiply( SysMatrix , divergence , miu , lamda , residual ) ;
%         beta = ( r_next' * r_next ) / ( r_previous' * r_previous ) ;  
%         d = r_next + beta * d ; 
%         
%         local_e = LocalError(img,img_previous);        
%         loss = norm( residual , 2) ;   % objective function
%         disp ( ['CG iteration: ', num2str(times),';   |    CG_local_e: ', num2str(local_e), ';   |    CG_Loss: ', num2str(loss)] ) ;
%         times = times + 1 ; 
%         img_previous = img ; 
%     end
% %     img = img + img_previous_iterative ;
% end
% 
% function output = matrixMultiply( SysMatrix , divergence , miu , lamda , colvect )
%     output = miu * (SysMatrix') * ( SysMatrix * colvect ) + lamda * divergence * colvect ;
% end


function img = cgls4TV ( SysMatrix, divergence, b_CG , iter_CG , miu , lamda, img_previous )
% using previous result to compute the search direction      
%     
%     img_previous = zeros(size(SysMatrix,2),1) ;     % initialization with previous result
%     b_CG = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous_iterative ) ; 
    
    % because the parameter matrix is miu*A'*A + lamda*(gx'*gx+gy'*gy), so
    % it is symetric and definite. To reduce the computation, here omit the
    % transpose operation on the parameter matrix. Since A is super large,
    % here i split the matrix, letting each component multiply separately.

    times = 1 ; 
    threshold = 1e-6 ;
    r0 = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img_previous ) ;   % initial residual for conjugated gradient algorithm ( CG ), which is negative to the gradient
    d = r0 ;    % initial search direction for conjugated gradient algorithm ( CG )     
    gama = norm(r0)^2 ;
    local_e = 100 ;   % initial
    resNE = 100 ;
    while ( times <= iter_CG && resNE >= threshold )

    % conjugated gradient algorithm ( CG )

        Ad = matrixMultiply( SysMatrix , divergence , miu , lamda , d ) ; 
% there is no difference using r_previous or using d to compute the alpha
%         alpha = ( d' * r_previous ) / ( d' * Ad + eps) ;       
             
        alpha = gama / ( d' * Ad + eps) ;
        img = img_previous + alpha * d ;
        residual = b_CG - matrixMultiply( SysMatrix , divergence , miu , lamda , img ) ; 
%         if ( mod ( times , 5 ) == 1 )
        r_next = residual ;          % there is no difference directly using residual or using iterative framework to compute the r_next
%         else 
%             r_next = r_previous - alpha * Ad ;      
%         end 
        gama1 = gama ;
        gama = norm(r_next)^2 ;
        beta = gama / ( gama1 + eps) ;  
        d = r_next + beta * d ; 
        
%         local_e = LocalError(img,img_previous);        
        gradient_loss = norm( residual ) ;   % objective function
%         loss = 0.5 * img' * matrixMultiply( SysMatrix , divergence ,  miu , lamda , img ) - b_CG' * img ; 
        resNE = norm(residual) / norm(r0) ;
%         disp ( ['CG iteration: ', num2str(times ),';   |    CG_local_e: ', num2str(local_e), ';   |    CG_resNE: ', num2str(resNE), ';   |    CG_Residual: ', num2str(gradient_loss), ';   |    CG_Loss: ', num2str(loss)] ) ;
        disp ( ['CG iteration: ', num2str(times ), ';   |    CG_resNE: ', num2str(resNE), ';   |    CG_Residual: ', num2str(gradient_loss)] ) ;
        img_previous = img ;
        
        times = times + 1 ;
    end
end

function output = matrixMultiply( SysMatrix , divergence ,  miu , lamda , colvect )
    output = miu * (SysMatrix') * ( SysMatrix * colvect ) + lamda * divergence * colvect ;
% output = miu * (SysMatrix') * ( SysMatrix * colvect ) ;
end



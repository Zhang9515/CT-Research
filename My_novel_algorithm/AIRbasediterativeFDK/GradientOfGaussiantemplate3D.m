function G = GradientOfGaussiantemplate3D ( ps )
    % 2019/01/02 by ZXZ
    % input: ps : size of gaussian patch; output: G : 3-ps-ps-ps GOG convolution template
    % 
    
    patchsize = [ps , ps , ps] ; 
    patchsizeG = [ 3 , patchsize ] ;
    
    G = zeros ( patchsizeG ) ; 
    Gx = zeros(patchsize);  
    Gy = zeros(patchsize);
    Gz = zeros(patchsize);
    sigma = 1 ; W = 0 ;
    for i = 1 : ps
        for j = 1 : ps
            for k = 1 : ps
                x = i - ps/2 - 0.5 ; y = j - ps/2 - 0.5; z = k - ps/2 - 0.5; 
                Gp = (-1/ ( (2*pi)^(3/2) *sigma^3)) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ; 
                Gx( i , j , k ) = Gp * (-x/sigma^2) ;  %
                Gy( i , j , k ) = Gp  * (-y/sigma^2);  %
                Gz( i , j , k ) = Gp  * (-z/sigma^2);  %
                W = W + Gp ;      % normalization factor
            end
        end
    end

    G(1 , : , : , : ) = Gx/W ;
    G(2 , : , : , : ) = Gy/W ;
    G(3 , : , : , : ) = Gz/W ;
    
end
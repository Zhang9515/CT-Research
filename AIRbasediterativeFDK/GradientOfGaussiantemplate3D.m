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
    sigma = 1 ; 
    for i = 1 : ps
        for j = 1 : ps
            for k = 1 : ps
                x = i - ( ps -1)/2 ; y = j - ( ps -1)/2 ; z = k - ( ps -1)/2 ; 
                Gx( i , j , k ) = (-1/ ( (2*pi)^(3/2) *sigma^3)) * (1-x/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Gy( i , j , k ) = (-1/ ( (2*pi)^(3/2) *sigma^3)) * (1-y/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Gz( i , j , k ) = (-1/ ( (2*pi)^(3/2) *sigma^3)) * (1-z/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
            end
        end
    end
    G(1 , : , : , : ) = Gx ;
    G(2 , : , : , : ) = Gy ;
    G(3 , : , : , : ) = Gz ;

end
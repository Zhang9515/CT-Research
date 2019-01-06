function H = HessianLOGtemplate3D( ps )
% 2018/11/12
    % input: size of patch( constant);  output: 3 - 3 - ps - ps - ps
    patchsize = [ps , ps , ps] ; 
    patchsizeH = [ 3 , 3 , patchsize ];
    
    H = zeros ( patchsizeH ) ; 
    Hx = zeros(patchsize);  
    Hy = zeros(patchsize);
    Hz = zeros(patchsize);
    Hxy = zeros(patchsize);
    Hxz = zeros(patchsize);
    Hyz = zeros(patchsize);
    sigma = 1 ; 
    for i = 1 : ps
        for j = 1 : ps
            for k = 1 : ps
                x = i - ( ps -1)/2 ; y = j - ( ps -1)/2 ; z = k - ( ps -1)/2 ; 
                Hx( i , j , k ) = (-1/ ( (2*pi)^3 *sigma^5)) * (1-x^2/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Hy( i , j , k ) = (-1/ ( (2*pi)^3 *sigma^5)) * (1-y^2/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Hz( i , j , k ) = (-1/ ( (2*pi)^3 *sigma^5)) * (1-z^2/sigma^2) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Hxy( i , j , k ) = ( x * y / ( (2*pi)^3 * sigma^7) ) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Hxz( i , j , k ) = ( x * z / ( (2*pi)^3 * sigma^7) ) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
                Hyz( i , j , k ) = ( y * z / ( (2*pi)^3 * sigma^7) ) * exp(-(x^2+y^2+z^2)/(2*sigma^2)) ;
            end
        end
    end
    H(1 , 1 , : , : , : ) = Hx ;
    H(1 , 2 , : , : , : ) = Hxy ;
    H(1 , 3 , : , : , : ) = Hxz ;
    H(2 , 1 , : , : , : ) = Hxy ;
    H(2 , 2 , : , : , : ) = Hy ;
    H(2 , 3 , : , : , : ) = Hyz;
    H(3 , 1 , : , : , : ) = Hxz ;
    H(3 , 2 , : , : , : ) = Hyz ;
    H(3 , 3 , : , : , : ) = Hz ;
end



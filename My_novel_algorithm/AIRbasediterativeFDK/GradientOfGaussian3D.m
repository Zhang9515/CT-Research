function gradient3D = GradientOfGaussian3D( image , ps , sizeofData )
    % 2019/01/02
    % assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
    % X : from left to right; Y : from down to up; Z : from near to far
    %input: image: N(=X-Y-Z)-1, ps: length of patch; output: gradient3D: N( =X*Y*Z )-3 
    image = reshape( image , sizeofData ) ; 
    gradient3D = zeros ( prod( sizeofData ) , 3 ) ;
    G = GradientOfGaussiantemplate3D( ps ) ;
    
    for i = 1 : 3
        gradient3D( : , i ) = reshape ( convn( image , squeeze(G( i , : , : , : )) , 'same' ) , prod( sizeofData ) , 1 ) ;             
    end

end
function gradient3D = GradientOfGaussian3D( image )
    % 2019/01/02
    % assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
    % X : from left to right; Y : from down to up; Z : from near to far
    %input: image: X-Y-Z, ps: length of patch; output: gradient3D: N( =X*Y*Z )-3 
    [ Lx, Ly, Lz ] = size ( image ) ;
    gradient3D = zeros ( Lx * Ly * Lz , 3 ) ;
    ps = 3 ;
    G = GradientOfGaussiantemplate3D( ps ) ;
    
    for i = 1 : 3
        gradient3D( : , i ) = reshape ( convn( image , G( i , : , : , : ) , 'same' ) , Lx * Ly * Lz , 1 ) ;             
    end

end
function gradient3D = Sobel3D( image, sizeofData )
    % input: image: N(X*Y*Z)*1, sizeofData: 1*3 vector; output: N(X*Y*Z) * 3
    image = reshape( image , sizeofData ) ; 
    gradient3D = zeros ( prod( sizeofData ) , 3 ) ;

    % sobel 3D operator
    sobelx = zeros( 3 , 3 , 3) ;
    sobelx(1,:,:) = [1 2 1; 2 12 2; 1 2 1] ; sobelx(3,:,:) = -sobelx(1,:,:) ;
    sobely = zeros( 3 , 3 , 3) ;
    sobely(:,1,:) = [1 2 1; 2 12 2; 1 2 1] ; sobely(:,3,:) = -sobelx(:,1,:) ;
    sobelz = zeros( 3 , 3 , 3) ;
    sobelz(:,1,:) = [1 2 1; 2 12 2; 1 2 1] ; sobelz(:,3,:) = -sobelx(:,1,:) ;
    
    gradient3D( : , 1 ) = reshape ( convn( image , sobelx , 'same' ) , prod( sizeofData ) , 1 ) ;          
    gradient3D( : , 2 ) = reshape ( convn( image , sobely , 'same' ) , prod( sizeofData ) , 1 ) ;
    gradient3D( : , 3 ) = reshape ( convn( image , sobelz , 'same' ) , prod( sizeofData ) , 1 ) ;
  
end
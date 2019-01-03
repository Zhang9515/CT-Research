function output = GOGtranspose3D( gradient, ps , sizeofData )
    % 2019/01/02
    % assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
    % X : from left to right; Y : from down to up; Z : from near to far
    % input: gradient3D: N(=X*Y*Z )-3 ; output:  N(=X*Y*Z )-1
    
    G = GradientOfGaussiantemplate3D( ps ) ;
    output = zeros ( sizeofData ) ;
    
     % actually we are computing the adjoint of hessian3D matrix
    G_transpose = G ;
    % exchange the corresponding row/column between the specific axis
        num = ceil (( ps - 1 ) / 2) ;
        c = G_transpose( : , 1 : num , : , : ) ; 
        G_transpose( : , 1 : num , : , : ) = G_transpose( : , ps - num +1 , : , : ) ; 
        G_transpose( : , ps - num +1 , : , : ) = c ;
        
        c = G_transpose( : , : , 1 : num , : ) ;
        G_transpose( : , : , 1 : num , : ) = G_transpose( : , : , ps - num +1 , : ) ; 
        G_transpose( : , : , ps - num +1 , : ) = c ; 
        
        c = G_transpose( : , : , : , 1 : num ) ;
        G_transpose( : , : , : , 1 : num ) = G_transpose( : , : , : , ps - num +1 ) ; 
        G_transpose( : , : , : , ps - num +1 ) = c ;
    
   for i =1 : 3
        output = output + convn( gradient( : , : , : , i ) , G_transpose ( ) , 'same' ) ;
   end
   
   output = reshape( output , prod( sizeofData ) ,1 ) ;
      
end
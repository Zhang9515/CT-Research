function result = HessianLOG3D( input )
% 2018/11/14 
% Input: x-y-z-dim vector; Output: 3-3-N(=x*y*z) matrix
    [ Lx , Ly, Lz ] = size( input ) ;
    result = zeros( 3 , 3 , Lx * Ly * Lz ) ;
    ps = 3 ;
    H = HessianLOGtemplate3D( ps ) ;
    
    for i = 1 : 3
        for j = 1 : 3
            result( i , j , : ) = reshape ( convn( input , H( i , j , : , : , : ) , 'same' ) , Lx * Ly * Lz , 1 ) ;
        end
    end
  
end












function Fai = ObjectiveFunction( ProjectionData, AU , U, OrderofHS, Tao)
% U is the matrix of image
    datalength = numel( U ) ;

    HSLOG3D = HessianLOG3D( U ) ; 
    for i = 1 : datalength
        svd_pixel = svd ( squeeze( HSLOG3D( : , : , i ) ) ) ;
        Hsnorm = Hsnorm + sum( svd_pixel.^OrderofHS )^( 1 / OrderofHS ) ;        
    end
    Fai = 0.5 * (ProjectionData - AU)' * (ProjectionData - AU) + Tao * Hsnorm ;

end
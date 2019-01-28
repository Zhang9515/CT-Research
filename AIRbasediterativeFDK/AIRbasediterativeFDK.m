%% 2019/01/02 by ZXZ
% AIR based iterative method














%% TV denoising ( for 3D mainly )

Input0 = phantom3d( 'Modified Shepp-Logan' , 128 ) ;
sizeofData = [ 128 , 128 , 128 ] ;
patchsize = 3 ; 
miu = 2 ;     % related to step size of each iteration
lamda = 1 ; 

X = zeros ( prod(sizeofData) , 1 ) ;
Z = zeros ( prod(sizeofData) , 3 ) ;
U = zeros ( prod(sizeofData) , 3 ) ;

for i = 1 : 128 
    Input( : , : , i ) = imnoise ( squeeze( Input0( : , : , i ) )  ,'salt & pepper', 0.02 ) ;
end
figure,imshow3Dfull ( Input , [ 0 1 ],'grey')
reference = reshape ( Input , prod(sizeofData) , 1 ) ;

iterationNumMax = 100 ; 
RMSE_MOD = zeros( iterationNumMax ,1 ) ;
ObjectValue = zeros( iterationNumMax ,1 ) ;

for iter = 1 : iterationNumMax
    disp(['iteration: ',num2str(iter)]) ;
    
    b = reference + miu * GOGtranspose3D( Z - U, patchsize , sizeofData ) ; 
    Residual = X + GOGtranspose3D ( GradientOfGaussian3D ( X , patchsize , sizeofData ) , patchsize , sizeofData ) - b ;
    alpha = dot ( Residual ,Residual ) / ( dot ( Residual , Residual + GOGtranspose3D ...
        ( GradientOfGaussian3D ( Residual , patchsize , sizeofData ) , patchsize , sizeofData ) ) + 1e-6 ) ;
    
    X = X - alpha * Residual ; 
    X = max ( X , 0 ) ;      % non-negative constrained
    X = min ( X , 1 ) ;      % non-negative constrained
    
    gradientofX = GradientOfGaussian3D ( X , patchsize , sizeofData ) ; 
    
    Z = isotropicISTA ( gradientofX + U , lamda / miu ) ; 
    
    U = U + gradientofX - Z ; 
%     RMSE_MOD(iter) = RMSE( X , reference ) ; 
%     disp(['RMSE: ',num2str(RMSE_MOD(iter))]) ;
    ObjectValue(iter) = 0.5 * norm( X - reference , 2 ) + lamda * norm( gradientofX , 1) ; 
    disp(['ObjectValue: ',num2str(ObjectValue(iter))]) ;
end

Display = reshape ( X , sizeofData ) ; 
figure,imshow3Dfull ( Display , [ 0 1 ],'grey')






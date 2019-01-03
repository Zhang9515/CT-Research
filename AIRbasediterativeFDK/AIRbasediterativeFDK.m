%% 2019/01/02 by ZXZ
% AIR based iterative method














%% TV denoising ( for 3D mainly )

Input0 = phantom3d( 'Modified Shepp-Logan' , 128 ) ;
sizeofData = [ 128 , 128 , 128 ] ;
patchsize = 3 ; 
miu = 1 ; 
lamda = 1 ; 

X = zeros ( prod(sizeofData) , 1 ) ;
Z = zeros ( prod(sizeofData) , 3 ) ;
U = zeros ( prod(sizeofData) , 3 ) ;

for i = 1 : 128 
    Input( : , : , i ) = imnoise ( squeeze( Input0( : , : , i ) )  ,'salt & pepper',0.02 ) ;
end
% figure,imshow3Dfull ( Input , [ 0 1 ])
reference = reshape ( Input , prod(sizeofData) , 1 ) ;

iterationNumMax = 20 ; 

for iter = 1 : iterationNumMax
    disp(['iteration:',num2str(iter)]) ;
    b = reference + miu * GOGtranspose3D( Z - U, patchsize , sizeofData ) ; 
    Residual = X + GOGtranspose3D ( GradientOfGaussian3D ( X , patchsize , sizeofData ) , patchsize , sizeofData ) - b ;
    alpha = dot ( Residual ,Residual ) / ( dot ( Residual , Residual + GOGtranspose3D ...
        ( GradientOfGaussian3D ( Residual , patchsize , sizeofData ) , patchsize , sizeofData ) ) + 1e-6 ) ;
    X = X - alpha * Residual ; 
    gradientofX = GradientOfGaussian3D ( X , patchsize , sizeofData ) ; 
    Z = isotropicISTA ( gradientofX + U , lamda / miu ) ; 
    U = U + gradientofX - Z ; 
end

Display = reshape ( X , sizeofData ) ; 
figure,imshow3Dfull ( Display , [ 0 1 ])






%% 2019/01/02 by ZXZ
% AIR based iterative method
tic
EPS = 1e-15 ; 













%% TV denoising ( for 3D mainly )

Input0 = phantom3d( 'Modified Shepp-Logan' , 128 ) ;
sizeofData = [ 128 , 128 , 128 ] ;
Input0flat = reshape(Input0 , prod(sizeofData) , 1 ) ; 
patchsize = 3 ; 
miu = 2 ;     % related to step size of each iteration
lamda = 0.5 ; 

X = zeros ( prod(sizeofData) , 1 ) ;       % corresponding to size of image
Z = zeros ( prod(sizeofData) , 3 ) ;       % corresponding to size of image gradient
U = zeros ( prod(sizeofData) , 3 ) ;       % corresponding to size of image gradient

for i = 1 : 128 
    Input( : , : , i ) = imnoise ( squeeze( Input0( : , : , i ) )  ,'salt & pepper',0.02 ) ;
end
Input = Input + ( Input - Input0 ) * 0.5 ;    % scale down the salt & pepper noise

% figure,imshow3Dfull ( Input , [ 0 1 ], 'grey')
reference = reshape ( Input , prod(sizeofData) , 1 ) ;

iterationNumMax = 5 ; 
RMSE_MOD = zeros( iterationNumMax ,1 ) ;
ObjectValue = zeros( iterationNumMax ,1 ) ;

for iter = 1 : iterationNumMax
    disp(['iteration: ',num2str(iter)]) ;
    
    b = reference + miu * GOGtranspose3D( Z - U, patchsize , sizeofData ) ; 
    Residual = X + GOGtranspose3D ( GradientOfGaussian3D ( X , patchsize , sizeofData ) , patchsize , sizeofData ) - b ;
    % alpha = rTr/rTpr 
    alpha = dot ( Residual ,Residual ) / ( dot ( Residual , Residual + miu * GOGtranspose3D ...
        ( GradientOfGaussian3D ( Residual , patchsize , sizeofData ) , patchsize , sizeofData ) ) + EPS ) ;
    
    X = X - alpha * Residual ; 
    X = max ( X , 0 ) ;      % non-negative constrained
    X = min ( X , 1 ) ;      % max-value constrained
    
    gradientofX = GradientOfGaussian3D ( X , patchsize , sizeofData ) ; 
    
    Z = isotropicISTA ( gradientofX + U , lamda / miu ) ; 
    
    U = U + gradientofX - Z ; 
    RMSE_MOD(iter) = RMSE( X , Input0flat ) ; 
    disp(['RMSE: ',num2str(RMSE_MOD(iter))]) ;
    ObjectValue(iter) = 0.5 * norm( X - reference , 2 ) + lamda * norm( gradientofX , 1) ; 
    disp(['ObjectValue: ',num2str(ObjectValue(iter))]) ;
end

Display = reshape ( X , sizeofData ) ; 
figure,imshow3Dfull ( Display , [ 0 1 ] , 'grey')
RMSEDisplay = RMSE3d(Input0, Display) ;
disp(['RMSEDisplay: ',num2str(RMSEDisplay)]) ;
toc




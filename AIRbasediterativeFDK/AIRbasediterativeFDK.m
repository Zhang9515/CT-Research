%% 2019/01/02 by ZXZ
% AIR based iterative method














%% TV denoising ( for 3D mainly )

Input0 = phantom3d( 'Modified Shepp-Logan' , 128 ) ;
sizeofData = [ 128 , 128 , 128 ] ;
patchsize = 3 ; 

for i = 1 : 128 
    Input( : , : , i ) = imnoise ( squeeze( Input0( : , : , i ) )  ,'salt & pepper',0.02 ) ;
end
% figure,imshow3Dfull ( Input , [ 0 1 ])
reference = Input ;

output = GradientOfGaussian3D ( Input , patchsize ) ;
output = reshape ( output( : , 1) , sizeofData ) ;
% figure,imshow3Dfull ( output , [  ])

% iterationNumMax = 20 ; 
% 
% for iter = 1 : iterationNumMax
%     
%     
%     
%     
%     
%     
%     
% end









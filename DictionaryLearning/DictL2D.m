%2019/03/21
% implement primary dictionary learning with OMP & K_SVD

load 'E:\ZXZ\Data\trial2D'
figure,imshow(trial2D,[0 0.5])
imgsize = size( trial2D ) ;
patchsize = [5 , 5] ; 
slidestep = [3 , 3] ;
patchset = ExtractPatch2D ( trial2D , patchsize , slidestep) ; 
Image2D = PatchSynthesis ( patchset , patchsize , slidestep, imgsize ) ;
figure,imshow(Image2D,[0 0.5])
% Alpha = omp( Dictionary , Xintm , Dictionary' * Dictionary , sparsity ) ;



 





















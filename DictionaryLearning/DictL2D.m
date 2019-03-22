%2019/03/21
% implement primary dictionary learning with OMP & K_SVD

load 'E:\ZXZ\Data\trial2D'
size = size( trial2D ) ;



Alpha = omp( Dictionary , Xintm , Dictionary' * Dictionary , sparsity ) ;

























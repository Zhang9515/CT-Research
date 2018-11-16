% 2018/11/12
% assume result U as column vector, which is tranformed from 3D volumn( x , y , z)
% X : from left to right; Y : from down to up; Z : from near to far
% 
clear all
tic
U = ones( 128, 128, 64 ) ;   % initialize U by result of FDK
[ Lx , Ly , Lz ] = size ( U ) ;
omega = ones ( 9 , Lx * Ly * Lz ) ;    % initalize omega0
ps = 3;
H = HessianLOGtemplate3D( ps ) ;

omega_stack = repmat ( omega , 1, 1, ps, ps, ps) ;

H_omega = zeros( Lx * Ly * Lz , 1 ) ; 

% actually we are computing the adjoint of hessian3D matrix
H_transpose = H ;
for num = 1 : ( ps - 1 ) / 2
    H_transpose( : , [num,3] , : , : ) = H_transpose( : , [num,1] , : , : ) ; 
    H_transpose( : , : , [num,3] , : ) = H_transpose( : , : , [num,1] , : ) ; 
    H_transpose( : , : , : , [num,3] ) = H_transpose( : , : , : , [num,1] ) ; 
end

for UJindex = 1 : Lx * Ly * Lz

            H_omega( UJindex ) = sum(sum(sum(sum ( squeeze( omega_stack( : , UJindex , : , : , : ) ) .* H_transpose ) ) ) ) ; 
              
end       

toc








































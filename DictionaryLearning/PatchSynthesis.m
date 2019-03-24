function Image2D = PatchSynthesis ( patchset , patchsize , slidestep, imgsize )
    % 2019/03/22 by ZXZ
    % synthsize patch into image
    % output of function, size is same as the size of image from which the patches are extracted 
    Image2D = zeros ( imgsize ) ;      
    % setup variable: Count to record the times of each pixels are used
    Count = zeros( imgsize ) ;
    % input: dictionary in which the patches are restored
    [ ms , ns ] = size( patchset ) ;
    % a all-ones template used when the patch cover the area
    Model = ones( ms , ns ) ;
    % the number of patches extending in two directions, thus the whole
    % number of patches is setm * setn
    setm = ceil(( imgsize(1) - patchsize(1) ) / slidestep(1) + 1) ;
    setn = ns / setm ;
    % cover the patches into the image domain. To consider the boundary
    % condition, the operations are carried out via 4 different situations.
    for i = 1 :  setm - 1
        for j = 1 : setn - 1
            Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset ( : , ( i - 1 ) * setn + j ) , patchsize ) ;
            Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( Model ( : , ( i - 1 ) * setn + j ) , patchsize ) ;
        end
    end
    
    for i = 1 :  setm - 1
        Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset ( : , i * setn ) , patchsize ) ;
        Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( Model ( : , i * setn ) , patchsize ) ;
    end
    
    for j = 1 : setn - 1
        Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset ( : , ( setm - 1 ) * setn + j ) , patchsize ) ;
        Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( Model ( : , ( setm - 1 ) * setn + j ) , patchsize ) ;
    end
    Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset ( : , setm * setn ) , patchsize ) ;
    Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( Model ( : , setm * setn ) , patchsize ) ;
    Image2D = Image2D ./ Count ;
end
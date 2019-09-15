function Image2D = PatchSynthesis ( patchset_RemoveDC , patchset , patchsize , slidestep, imgsize, alternate )
    % 2019/03/22 by ZXZ
    % image represents the original image before the process
    % synthsize patch into image
    defaults = {'AddDC','NoAddDC'} ;
    idx = find(strncmp(defaults, alternate, 5 ) ) ;
    
    % output of function, size is same as the size of image from which the patches are extracted 
    Image2D = zeros ( imgsize ) ;      
    % setup variable: Count to record the times of each pixels are used
    Count = zeros( imgsize ) ;
    % input: dictionary in which the patches are restored
    [ ms , ns ] = size( patchset_RemoveDC ) ;
    % a all-ones template used when the patch cover the area
    Model = ones( patchsize ) ;
    % the number of patches extending in two directions, thus the whole
    % number of patches is setm * setn
    setm = ceil(( imgsize(1) - patchsize(1) ) / slidestep(1) + 1) ;
    setn = ns / setm ;
    % cover the patches into the image domain. To consider the boundary
    % condition, the operations are carried out via 4 different situations.
    if (idx ==1)
            for i = 1 :  setm - 1
                for j = 1 : setn - 1
                    Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset_RemoveDC ( : , ( i - 1 ) * setn + j ) , patchsize ) +... 
                    mean2 (patchset( : , ( i - 1 ) * setn + j ) ) ;
                    Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = Model ;
                end
            end

            for i = 1 :  setm - 1
                Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset_RemoveDC ( : , i * setn ) , patchsize ) +... 
                    mean2 (patchset( : , i * setn ) ) ;
                Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = Model ;
            end

            for j = 1 : setn - 1
                Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset_RemoveDC ( : , ( setm - 1 ) * setn + j ) , patchsize ) +... 
                    mean2 (patchset( : , ( setm - 1 ) * setn + j ) ) ;
                Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = Model ;
            end
            Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset_RemoveDC ( : , setm * setn ) , patchsize ) +... 
                    mean2 (patchset( : , setm * setn ) ) ;
            Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = Model ;
    else
            for i = 1 :  setm - 1
                for j = 1 : setn - 1
                    Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset_RemoveDC ( : , ( i - 1 ) * setn + j ) , patchsize );
                    Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = Model ;
                end
            end

            for i = 1 :  setm - 1
                Image2D ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset_RemoveDC ( : , i * setn ) , patchsize );
                Count ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = Model ;
            end

            for j = 1 : setn - 1
                Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = reshape( patchset_RemoveDC ( : , ( setm - 1 ) * setn + j ) , patchsize );
                Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) = Model ;
            end
            Image2D ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = reshape( patchset_RemoveDC ( : , setm * setn ) , patchsize );
            Count ( imgsize(1) - patchsize(1) + 1  : imgsize(1) , imgsize(2) - patchsize(2) + 1 : imgsize(2) ) = Model ;
    end
    
    Image2D = Image2D ./ Count ;
end
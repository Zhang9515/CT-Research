function patchset = ExtractPatch2D ( image , patchsize , slidestep, alternate )
    % 2019/03/22 by ZXZ
    % extract patch from image
    defaults = {'RemoveDC','NoRemoveDC'} ;
    idx = find(strncmp(defaults, alternate, 8 ) ) ;
    
    [ m , n ] = size( image ) ;
    setm = ceil(( m - patchsize(1) ) / slidestep(1) + 1) ;
    setn = ceil(( n - patchsize(2) ) / slidestep(2) + 1 ) ;
    patchset = zeros ( prod(patchsize) , setm * setn ) ;
    % operate DC remove within each patches before being extrated
    if (idx ==1)
            for i = 1 :  setm - 1
                for j = 1 : setn - 1
                    patchset ( : , ( i - 1 ) * setn + j ) = reshape( removeDC ( image ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) ) , prod(patchsize)  , 1 ) ;     
                end
            end

            for i = 1 :  setm - 1
                patchset ( : , i * setn ) = reshape( removeDC ( image ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , n - patchsize(2) + 1 : n ) ) , prod(patchsize)  , 1 ) ;     
            end

            for j = 1 : setn - 1
                patchset ( : , ( setm - 1 ) * setn + j ) = reshape( removeDC ( image ( m - patchsize(1) + 1  : m , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) ) , prod(patchsize)  , 1 ) ; 
            end

            patchset ( : , setm * setn ) = reshape( removeDC ( image ( m - patchsize(1) + 1  : m , n - patchsize(2) + 1 : n ) ) , prod(patchsize)  , 1 ) ;
    % no operation of DC removal
    elseif(idx ==2)
            for i = 1 :  setm - 1
                for j = 1 : setn - 1
                    patchset ( : , ( i - 1 ) * setn + j ) = reshape( image ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) ) , prod(patchsize)  , 1 ) ;     
                end
            end

            for i = 1 :  setm - 1
                patchset ( : , i * setn ) = reshape( image ( 1 + ( i - 1 ) * slidestep(1) : patchsize(1) + ( i - 1 ) * slidestep(1) , n - patchsize(2) + 1 : n )  , prod(patchsize)  , 1 ) ;     
            end

            for j = 1 : setn - 1
                patchset ( : , ( setm - 1 ) * setn + j ) = reshape( image ( m - patchsize(1) + 1  : m , 1 + ( j - 1 ) * slidestep(2) : patchsize(2) + ( j - 1 ) * slidestep(2) )  , prod(patchsize)  , 1 ) ; 
            end

            patchset ( : , setm * setn ) = reshape( image ( m - patchsize(1) + 1  : m , n - patchsize(2) + 1 : n ) , prod(patchsize)  , 1 ) ;
    
     end
    

end
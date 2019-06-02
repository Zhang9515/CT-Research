function patchset = AdaptiveExtractPatch2D ( image , Patchsize , alternate )
    % 2019/05/31 by ZXZ
    % extract patch from image adaptively
    % input: Patchsize is the patch size map
    % patchset is cell tuple of patch with various size
    defaults = {'RemoveDC','NoRemoveDC'} ;
    idx = find(strncmp(defaults, alternate, 8 ) ) ;
    
    [ m , n ] = size( image ) ;
    patchset = cell( m*n , 1 ) ;
    % operate DC remove within each patches before being extrated
    if (idx ==1)
            for index_x = 1 : n
                for index_y = 1 : m
                    patchsize = Patchsize ( index_y , index_x ) ;
                    if ( patchsize)
                end
            end          
    % no operation of DC removal
    elseif(idx ==2)
           
    end
    
end
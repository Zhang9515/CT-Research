function patchset = AdaptiveExtractPatch2D ( image , Patchsize , alternate )
    % 2019/06/02 by ZXZ
    % extract patch from image adaptively
    % input: Patchsize is the patch size map transformed into a column
    % vector
    % patchset is cell tuple of patch with various size
    defaults = {'RemoveDC','NoRemoveDC'} ;
    idx = find(strncmp(defaults, alternate, 8 ) ) ;
    
    [ height , width] = size( image ) ;
    patchset = cell( height*width , 1 ) ;
    % operate DC remove within each patches before being extrated
    if (idx ==1)
            for index_x = 1 : width
                for index_y = 1 : height
                    Img_index_y = height - index_y + 1 ;
                    index_pix = index_x + (index_y-1) * width ;
                    patchsize = Patchsize ( index_pix ) ;
                    if ( patchsize ~= 0 )
                        shift = (patchsize-1) / 2 ;
                        patchset { index_pix } = reshape ( removeDC( image (Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) ) , patchsize*patchsize , 1) ;
                    end
                end
            end          
    % no operation of DC removal
    elseif(idx ==2)
           for index_x = 1 : width
                for index_y = 1 : height
                    Img_index_y = height - index_y + 1 ;
                    index_pix = index_x + (index_y-1) * width ;
                    patchsize = Patchsize ( index_pix ) ;
                    if ( patchsize ~= 0 )
                        shift = (patchsize-1) / 2 ;
                        patchset { index_pix } = reshape ( image (Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) , patchsize*patchsize , 1) ;
                    end
                end
            end   
    end  
end
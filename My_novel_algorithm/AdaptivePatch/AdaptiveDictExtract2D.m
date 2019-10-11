function versatileDictionary = AdaptiveDictExtract2D ( Ref_image , Patchsize , slidestep , patch_level , interval )
    % 2019/10/10 by ZXZ
    % extract patch from Ref_image adaptively
    % input: Patchsize is the patch size map transformed into a column
    % vector
    % patchset is cell tuple of patch with various size
    versatileDictionary = cell( patch_level ,1 ) ;
    
    [ height , width] = size( Ref_image ) ;
    setm = floor(( height - 1 ) / slidestep(1) + 1) ;
    setn = floor(( width - 1 ) / slidestep(2) + 1 ) ;
    
    for index_x = 1 : setn
        for index_y = 1 : setm
            Img_index_y = height - (index_y-1) * slidestep(1) ;
            Img_index_x = (index_x-1) * slidestep(2) +1;
            index_pix = Img_index_x + (Img_index_y-1) * width ;
            patchsize = Patchsize ( index_pix ) ;
            shift = (patchsize-1) / 2 ;
            if ( Img_index_y - shift <1)
                py1 = 1 ;
                py2 = 2 * shift +1 ;
            elseif ( Img_index_y + shift >height )
                py1 = height - 2 * shift ;
                py2 = height ;
            else
                py1 = Img_index_y - shift ;
                py2 = Img_index_y + shift ;
            end

            if (Img_index_x - shift <1)
                px1 = 1 ;
                px2 = 2 * shift +1 ;
            elseif ( Img_index_x + shift >width )
                px1 = width - 2 * shift ;
                px2 = width ;
            else
                px1 = Img_index_x - shift ;
                px2 = Img_index_x + shift ;
            end
            
            level = floor((patchsize-1)/(2 * interval) )+1 ;
            versatileDictionary { level }(:,end+1) = reshape ( Ref_image (py1 : py2 , px1 : px2 ) , patchsize*patchsize , 1) ;                 

        end
    end     
    for index = 1 : patch_level
        versatileDictionary{ index } = col_normalization( versatileDictionary{ index } ) ;
    end
end
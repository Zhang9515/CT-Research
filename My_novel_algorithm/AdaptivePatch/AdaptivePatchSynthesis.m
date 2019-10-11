function Image2D = AdaptivePatchSynthesis ( patchset_RemoveDC_tuple , patchset_tuple, imgsize, slidestep, alternate )
    % 2019/05/31 by ZXZ
    % image represents the original image before the process
    % synthsize patch into image
    defaults = {'AddDC','NoAddDC'} ;
    idx = find(strncmp(defaults, alternate, 5 ) ) ;
    
    % output of function, size is same as the size of image from which the patches are extracted 
    Image2D = zeros ( imgsize ) ;      
    % setup variable: Count to record the times of each pixels are used
    Count = zeros( imgsize ) ;
    % input: all patches corresponding to the pixels has been stored in the
    % tuple, some pixels on the boundry hasn't patches which is empty in
    % the corresponding location
%     if (idx ==1)
%             for index_y = 1: imgsize(1)
%                 for index_x = 1: imgsize(2)     
%                     Img_index_y = imgsize(1) - index_y + 1 ;
%                     pixel_index = index_x + ( index_y -1 ) * imgsize(2) ;
%                     patch_RemoveDC = patchset_RemoveDC_tuple{pixel_index} ;
%                     patch = patchset_tuple {pixel_index} ;
%                     if ( ~isempty( patch_RemoveDC ) )
%                         patchsize = sqrt( numel( patch_RemoveDC ) ) ;
%                         shift = (patchsize-1) / 2 ;
%                         Model = ones( patchsize , patchsize ) ;       % a all-ones template used when the patch cover the area
%                         Image2D ( Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) = reshape( patch_RemoveDC , patchsize , patchsize ) + mean2( patch ) ;
%                         Count ( Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) = Model ;
%                     end
%                 end
%             end
%     else
%             for index_y = 1: imgsize(1)
%                 for index_x = 1: imgsize(2)
%                     Img_index_y = imgsize(1) - index_y + 1 ;
%                     pixel_index = index_x + ( index_y -1 ) * imgsize(2) ;
%                     patch_RemoveDC = patchset_RemoveDC_tuple{pixel_index} ;
%                     if ( ~isempty( patch_RemoveDC ) )
%                         patchsize = sqrt( numel( patch_RemoveDC ) ) ;
%                         shift = (patchsize-1) / 2 ;
%                         Model = ones( patchsize , patchsize ) ;       % a all-ones template used when the patch cover the area
%                         Image2D ( Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) = reshape( patch_RemoveDC , patchsize , patchsize ) ;
%                         Count ( Img_index_y - shift : Img_index_y + shift , index_x - shift : index_x + shift ) = Model ;
%                     end
%                 end
%             end
%     end
    
    setm = floor(( imgsize(1) - 1 ) / slidestep(1) + 1) ;
    setn = floor(( imgsize(2) - 1 ) / slidestep(2) + 1 ) ;
    
    for index_x = 1 : setn
        for index_y = 1 : setm
            Img_index_y = imgsize(1) - (index_y-1) * slidestep(1) - 1 ;
            Img_index_x = (index_x-1) * slidestep(2) +1;
            index_pix = Img_index_x + (Img_index_y-1) * imgsize(2) ;
            patch_RemoveDC = patchset_RemoveDC_tuple{index_pix} ;
            patch = patchset_tuple {index_pix} ;
            patchsize = sqrt( numel( patch_RemoveDC ) ) ;
            shift = (patchsize-1) / 2 ;
            Model = ones( patchsize , patchsize ) ;       % a all-ones template used when the patch cover the area
            
            if ( Img_index_y - shift <1)
                py1 = 1 ;
                py2 = 2 * shift +1 ;
            elseif ( Img_index_y + shift >imgsize(1) )
                py1 = imgsize(1) - 2 * shift ;
                py2 = imgsize(1) ;
            else
                py1 = Img_index_y - shift ;
                py2 = Img_index_y + shift ;
            end

            if (Img_index_x - shift <1)
                px1 = 1 ;
                px2 = 2 * shift +1 ;
            elseif ( Img_index_x + shift >imgsize(2) )
                px1 = imgsize(2) - 2 * shift ;
                px2 = imgsize(2) ;
            else
                px1 = Img_index_x - shift ;
                px2 = Img_index_x + shift ;
            end
            
            if ( idx == 1) % operate DC remove within each patches before being extrated               
                Image2D (py1 : py2 , px1 : px2 ) = reshape( patch_RemoveDC , patchsize , patchsize ) + mean2( patch ) ;
            elseif ( idx == 2)  % no operation of DC removal
                Image2D ( py1 : py2 , px1 : px2 ) = reshape( patch_RemoveDC , patchsize , patchsize ) ;       
            end
            Count ( py1 : py2 , px1 : px2 ) = Model ;
        end
    end   
    Image2D = Image2D ./ Count ;
end
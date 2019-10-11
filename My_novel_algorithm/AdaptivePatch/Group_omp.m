function patchset_tuple = Group_omp( Dictionary , Xintm, patchsize, imgsize,  sparsity, patch_level , interval )
    patchset_tuple = cell( prod(imgsize) ,1) ;
    maxindex = patch_level ;
    containerDictionary = cell( maxindex,1) ;
    containerXin = cell( maxindex,1) ;
    % acquire the index of input and store them in the various containers corresponding
    % to their size
    for index = 1 : prod(imgsize)
        if (patchsize(index)~=0)
            conD = floor((patchsize(index)-1)/(2 * interval) )+1 ;
            containerDictionary{conD}(end+1) = index ;     % store the source pixel index of the patch 
        end
    end
    % acquire the column patch vector depending on the index and do the OMP in groups 
    for index_D = 1 : maxindex
        disp(['index_D: ', num2str(index_D)])
        dic = Dictionary{index_D} ;
        Xinsize = length(containerDictionary{index_D}) ;
        if (Xinsize~=0)
            ps = (1+(index_D-1)*interval)*2+1 ;
            containerXin{index_D} = zeros( ps*ps, Xinsize) ;
            for index_xin = 1 : Xinsize
                if ( ~isempty( Xintm{ containerDictionary{index_D}(index_xin) } ) )
                    containerXin{index_D}(:,index_xin) = Xintm{ containerDictionary{index_D}(index_xin)} ;
                end
            end
            Alpha = omp( dic , containerXin{index_D} , dic' * dic , sparsity ) ; 
            containerXin{index_D} = dic * Alpha ;
        end     
    end
    % summarize the results distributed in various containers into output
    % tuple
    for index_D = 1 : maxindex
        Xinsize = length(containerDictionary{index_D}) ;
        if (Xinsize~=0)
                for index_xin = 1 : Xinsize    
                    patchset_tuple{containerDictionary{index_D}(index_xin)} = containerXin{index_D}(:,index_xin) ;
                end
        end        
    end
end
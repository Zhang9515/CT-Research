%19/05/23 by ZXZ
% adaptively select the proper patch size for every pixels
% reference to a paper which introduced a method based on the gradient angle
% histogram
clear all
tic
pic = phantom(256) ;
[height , width] = size(pic) ;
picvector = Img2vec_Mat2Cpp2D( pic ) ;
GradientxMatrix = gradient2Dmatrix_x( height , width ) ;
GradientyMatrix = gradient2Dmatrix_y( height , width ) ;
GX = GradientxMatrix * double(picvector) ;
GY = GradientyMatrix * double(picvector) ;
GradAng = rad2deg(atan2( GY , GX )) ;
PatchSize = zeros( numel(pic), 1 ) ;       % store the patch size corresponding to each pixel
edges = -180 : 10 : 180 ;
threshold = 0 ;      % KL threshold
hist = zeros( numel(pic), numel(edges)-1 ) ;     % store the histogram of gradient angle
for pz = 1 : 5
    %
   
    %
    for x = 1+pz : width-pz
        for y = 1+pz : height-pz          % ordinary area
            %
            x = 42 ; y = 126 ;
            %
            center_index =  x + (y-1) * width ;
            if (pz == 1)
                index = kron(-pz:pz,ones(2*pz+1,1)) + kron(ones(1,2*pz+1),width * (-pz:pz)') + center_index;
                hist(center_index,:) = histcounts( GradAng ( index ) , edges ) ;
                Pprevious = hist(center_index,:)   / sum(hist(center_index,:) ) ;
                PatchSize( center_index ) = 2*pz + 1 ;
                KLprevious = 0 ;
            else
                index = [ (-pz : pz) + pz * width, (-pz : pz) - pz * width, -pz + (-pz:pz)*width, pz + (-pz:pz)*width] + center_index ;
                hist(center_index,:)  = hist(center_index,:)  + histcounts( GradAng ( index ) , edges ) ;
                Pcurrent = hist(center_index,:) / sum(hist(center_index,:) )  ;
                KLcurrent = kldiv( 1:(numel(edges)-1) , Pprevious + eps, Pcurrent + eps, 'sym') ;      % compute KL divergence between current and previous GradAng distribution
                Pprevious = Pcurrent ;
                if ( (KLcurrent - KLprevious) > threshold)     % if current KL divergence is larger enough than the previous, we select the current patch size
                    KLprevious = KLcurrent ;
                    PatchSize( center_index ) = 2*pz + 1 ;
                end
            end    
        end      
    end
    % area need special operation on the boundry, there are 8 cases
    % altogether, denoted by left, right, up, down
    % but actually what I do above has perfectly solve the boundry area,
    % and all pixels are covered by at least one patch.  
end
PatchSize = Vec2img_Cpp2Mat2D( PatchSize , height , width ) ;
figure,imshow(PatchSize,[])
toc

%2018/06/25 by ZXZ
% K-SVD for dictionary learning & OMP for sparse encoding
tic
%%  intialize the dictionary
pic = phantom(512);       % original figure
Display = zeros(sizepic);    % reconstructed figure

sizepic = size(pic);
picflat = pic(:);
patchsize=[16,16];
Cdct = gen_dct(prod(patchsize));      % Discret Cosine Transform
%Wdwt = gen_wave(length(picflat),8);    % Discret wavelet Transform( daubechies D8)

Dictionary = [eye(prod(patchsize)), Cdct];
% y=displayPatch(Dictionary);
% y = Dictionary * picflat;
% ind = find(abs(y)<0.05);           % discard part of data
% comRate = length(ind) / length(picflat);
% y(ind) =0;
% picflatrev = Dictionary' * y;
% picrev = reshape(picflatrev,64,64);
% figure,imshow(picrev,[])
% title(comRate)

%% patch
stride = 8;
Train = zeros(prod(patchsize),prod( floor((sizepic-patchsize)/stride)+1));
for i = 1: floor((sizepic(1)-patchsize(1))/stride)+1
    for j = 1: floor((sizepic(2)-patchsize(2))/stride)+1
        patch = pic(1 + (i -1)*stride : 1 + (i -1)*stride+patchsize(1)-1, 1 + (j -1)*stride : 1 + (j -1)*stride+patchsize(2)-1);
        Train(:,(i-1)*(floor((sizepic(2)-patchsize(2))/stride)+1)+j) = patch(:)';
    end
end    
% y=displayPatch(Train);
Dictionary = randn(size(Train,1),256);
Dictionary = Dictionary./ (repmat(sqrt(sum(Dictionary.^2)),size(Dictionary,1),1)+1e-10);

noIt = 50;
mssims = zeros( noIt ,1);
rmses = zeros( noIt ,1);
for it = 1:noIt
    Display = zeros(sizepic);
    disp(strcat('itreation:', 32, num2str(it)))
    %% L = 5 Cholesky-factorization ( norm of atom is 1)
    Lthreshold = 20 ;
    Represent = zeros( size(Dictionary,2) ,size(Train,2));
    Alpha0 = Dictionary' * Train;
    for patchindex = 1: size(Train,2)   
        residual = Train(:,patchindex);
        Index =[];
        L = [1];
        for Ls = 1:Lthreshold
            err = abs(Dictionary' * residual) ;
            err(Index) = -1;
            [~,Mindex]=max(err);

            if (Ls>1)
                w = L^-1 * Dictionary( :, Index)'* Dictionary( :, Mindex) ;
                L = [ L zeros( size(L,1),1); w' sqrt(1-w'*w)];          
            end

            Index = [Index Mindex];
            Alpha = Alpha0( Index ,patchindex);
            gama = Cholesky(L, Alpha );
            residual = Train(:,patchindex) - Dictionary( :, Index) * gama ;    
        end
        Represent( Index ,patchindex) = gama;
    end

    data = Dictionary * Represent;
    % y=displayPatch(data);
    %% standard SVD
    % for it = 1:noIt
    %     W = sparseapprox(X, D, 'mexOMP', 'tnz',s);
    %     R = X - D*W;
    %     for k=1:K   % 
    %         I = find(W(k,:));
    %         Ri = R(:,I) + D(:,k)*W(k,I);
    %         [U,S,V] = svds(Ri,1,'L');
    %         % U is normalized
    %         D(:,k) = U;
    %         W(k,I) = S*V';
    %         R(:,I) = Ri - D(:,k)*W(k,I);
    %     end    
    % end
    % 
    %% approximate SVD


        %W = sparseapprox(X, D, 'mexOMP', 'tnz',s);      
        ERR = Train - Dictionary*Represent;
        for k=1:size(Dictionary,2)
            I = find(Represent(k,:));
            Ri = ERR(:,I) + Dictionary(:,k)*Represent(k,I);
            dk = Ri * Represent(k,I)';
            dk = dk/(sqrt(dk'*dk) + 1e-10);  % normalize
            Dictionary(:,k) = dk;
            Represent(k,I) = dk'*Ri;
            ERR(:,I) = Ri - Dictionary(:,k)*Represent(k,I);
        end
        
        for n = 1:size(Train,2)
            row = floor( (n -1) / (floor((sizepic(2)-patchsize(2))/stride)+1) ) +1;
            column = n - (row -1) * (floor((sizepic(2)-patchsize(2))/stride)+1);
            Display( 1 + (row -1)* stride: 1 + (row -1)* stride + patchsize(1) -1, 1 + (column -1)* stride : 1 + (column -1)* stride + patchsize(1) -1) ...
                = Display( 1 + (row -1)* stride: 1 + (row -1)* stride + patchsize(1) -1, 1 + (column -1)* stride : 1 + (column -1)* stride + patchsize(1) -1) ...
                + reshape(data(:,n),patchsize);
        end

        overlapstore =zeros(sizepic);
        for i = 1: sizepic(1)
            for j = 1: sizepic(2)
        %         overlap = max (( min( patchsize(1) -1, i -1) + min ( patchsize(1) -1, sizepic(1) - i) +1 - patchsize(1) + 1), 0 ) * ...
        %             max (( min( patchsize(2) -1, j -1) + min ( patchsize(2) -1, sizepic(2) - j) +1 - patchsize(2) + 1), 0 );
                overlap_rowindex = max( ceil( ( i -patchsize(1) )/stride )+1, 1 ) : min( floor ( (i-1)/stride )+1, floor( ( sizepic(1)-patchsize(1))/stride ) +1 );
                overlap_columnindex = max( ceil( ( j -patchsize(2) )/stride )+1, 1 ) : min( floor ( (j-1)/stride )+1, floor( ( sizepic(2)-patchsize(2))/stride ) +1 );
                overlap = length(overlap_rowindex) * length(overlap_columnindex);
                overlapstore(i,j) = overlap;
                %disp(overlap);
                if (overlap ~=0)
                    Display(i , j) = Display(i , j)/overlap;
                end
            end
        end    
        [mssims(it), ~] = SSIM(pic, Display);    
        rmses(it) = RMSE(pic, Display);    
        figure(1)   % SSIM
        plot ( 1 : it , mssims ( 1  : it ) ) ;
        ylim ( [ 0 , 1 ] ) ;
        title( strcat('SSIM', 32, 'number of atoms:', 32, num2str(Lthreshold), 32, 'Stride', 32, num2str(stride), 32, 'Patch', 32, num2str(patchsize(1))))
        drawnow ; 
        figure(2)   % RMSE
        plot ( 1 : it , rmses ( 1  : it ) ) ;
        ylim ( [ 0 , ( 10 * rmses ( it ) ) ] ) ;
        title( strcat('RMSE', 32, 'number of atoms:', 32, num2str(Lthreshold), 32, 'Stride', 32, num2str(stride), 32, 'Patch', 32, num2str(patchsize(1))))
        drawnow;
end

%% restore


% for n = 1:size(Train,2)
%     row = floor( (n -1) / (floor((sizepic(2)-patchsize(2))/stride)+1) ) +1;
%     column = n - (row -1) * (floor((sizepic(2)-patchsize(2))/stride)+1);
%     Display( 1 + (row -1)* stride: 1 + (row -1)* stride + patchsize(1) -1, 1 + (column -1)* stride : 1 + (column -1)* stride + patchsize(1) -1) ...
%         = Display( 1 + (row -1)* stride: 1 + (row -1)* stride + patchsize(1) -1, 1 + (column -1)* stride : 1 + (column -1)* stride + patchsize(1) -1) ...
%         + reshape(data(:,n),patchsize);
% end
% 
% %overlapstore =zeros(sizepic);
% for i = 1: sizepic(1)
%     for j = 1: sizepic(2)
% %         overlap = max (( min( patchsize(1) -1, i -1) + min ( patchsize(1) -1, sizepic(1) - i) +1 - patchsize(1) + 1), 0 ) * ...
% %             max (( min( patchsize(2) -1, j -1) + min ( patchsize(2) -1, sizepic(2) - j) +1 - patchsize(2) + 1), 0 );
%         overlap_rowindex = max( ceil( ( i -patchsize(1) )/stride )+1, 1 ) : min( floor ( (i-1)/stride )+1, floor( ( sizepic(1)-patchsize(1))/stride ) +1 );
%         overlap_columnindex = max( ceil( ( j -patchsize(2) )/stride )+1, 1 ) : min( floor ( (j-1)/stride )+1, floor( ( sizepic(2)-patchsize(2))/stride ) +1 );
%         overlap = length(overlap_rowindex) * length(overlap_columnindex);
%         %overlapstore(i,j) = overlap;
%         %disp(overlap);
%         if (overlap ~=0)
%             Display(i , j) = Display(i , j)/overlap;
%         end
%     end
% end    
figure,imshow(Display,[0 1])
title( strcat('number of atoms:', 32, num2str(Lthreshold), 32, 'Stride', 32, num2str(stride), 32, 'Patch', 32, num2str(patchsize(1))))
%figure,imshow(overlapstore,[])
%y=displayPatch(Dictionary);


 toc
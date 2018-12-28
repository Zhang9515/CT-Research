%% 2018/12/28 by ZXZ
% repair FDK image in Fourier domain by using prior image

%% 3D-FFT
img = Diskphantom ( 129 ) ;
% figure, imshow3Dfull( img , [0 1] )

% imgF = fft(fft(fft(img, [], 1), [], 2), [], 3) ; 
imgF = fftshift ( fftn ( img ) ) ;
MaxImgF = max(max(max( imgF ))) ; 

imgFlog = 1 / log( 10e-6 +  abs( MaxImgF ) ) * log( 10e-6 +  abs( imgF ) ) ;  
figure, imshow3Dfull( imgFlog , [] )

%%  low-pass filter
[tlength , slength , zlength] = size( img ) ; 
for tindex = 1 : tlength
    for sindex = 1 : slength
        for zindex = 1 : zlength
                distance = sqrt ( ( tindex - 0.5 - tlength / 2 )^2 + ( sindex - 0.5 - slength / 2 )^2 ...
                    + ( zindex - 0.5 - zlength / 2 )^2) ; 
                if distance > 30
                    imgF( tindex , sindex , zindex ) = 0 ; 
                end
        end
    end
end

imgreconstruction = real ( ifftn ( ifftshift ( imgF ) ) ) ; 
figure, imshow3Dfull( imgreconstruction , [] )
%% 2D-FFT
% img = phantom ( 129 ) ;
% % figure, imshow( img , [0 1] )
% 
% imgf = fft2( img ) ;
% imgf = fftshift( imgf ) ;    % move the low-frequency part to the center of the frequency profile
% figure, imshow( imgf , [ ] )
% 
% imgfscale = imgf * 100 ; 
% Maximgf = abs( max(max( imgfscale )) ) ;
% imgflog = 1 / log( 10e-6 + Maximgf ) * log( 10e-6 + abs( imgfscale ) ) ;
% figure, imshow( imgflog , [ ] )

%%


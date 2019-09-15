%reconstructed image
m = fft2(Display(32,:,:));
k = fftshift(m);
I = log(abs(k));
kr = real(k);
ki = imag(k);
A = sqrt(kr.^2 + ki.^2);
% A = (A-min(min(A)))/(max(max(A)) - min(min(A)))*255;
figure,imshow(squeeze(I),[])
% original image
m1 = fft2(double(pic(32,:,:)));
k1 = fftshift(m1);
I1= log(abs(k1));
kr = real(k);
ki = imag(k);
A = sqrt(kr.^2 + ki.^2);
% A = (A-min(min(A)))/(max(max(A)) - min(min(A)))*255;
figure,imshow(squeeze(I1),[])
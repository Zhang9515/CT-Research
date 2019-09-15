function y=displayPatch(w,window)
[n,Patchnum]=size(w);
psize=sqrt(n);
Patline=ceil(sqrt(Patchnum));
y=ones(Patline*(psize+1),Patline*(psize+1));
for i=1:Patchnum
    row=floor((i-1)/Patline)*(psize+1);
    column=mod((i-1),Patline)*(psize+1);
    z=transpose(reshape(w(:,i),psize,psize));
   %z=(z-min(z(:)))/(max(z(:))-min(z(:)));
    y(row+1:row+psize,column+1:column+psize)=z;
end
figure;imshow(y,window);
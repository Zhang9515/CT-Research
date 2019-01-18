% 2019-01-08 by ZXZ
% nifti data, interplotion through height
clear all
% nii = load_nii( 'G:\CTcode\Data\lung cancer\CT190108.nii', [], 1) ;
nii = load_nii( 'G:\CTcode\Data\CTThoraxAdonomen.nii', [], 1) ;
origin_img = double ( nii.img ) ;
origin_img = ( origin_img + 1024 ) / 4096 ;

horizontal_resolution = nii.hdr.dime.pixdim ( 2 ) ;
vertical_resolution = nii.hdr.dime.pixdim ( 4 ) ;
origin_size = nii.hdr.dime.dim ( 2 : 4 ) ;
augmentRate = floor( vertical_resolution/horizontal_resolution ) ;

current_size = origin_size ; current_size ( 3 ) = ( origin_size(3) - 1 ) * augmentRate + 1 ; 
xx = horizontal_resolution + horizontal_resolution * ( 0 : (origin_size(1)-1) ) ;
yy = xx ; 
zz = vertical_resolution + vertical_resolution * ( 0 : (origin_size(3)-1) ) ; 
[ XX , YY , ZZ ] = meshgrid ( xx , yy , zz ) ;

zzprime = vertical_resolution + vertical_resolution/augmentRate * ( 0 : (current_size(3)-1) ) ; 
[ XXprime , YYprime , ZZprime ] = meshgrid ( xx , yy , zzprime ) ;

current_img = interp3( XX , YY , ZZ , origin_img , XXprime , YYprime , ZZprime , '*cubic' ) ; 
figure,imshow3Dfull( current_img , [0 1],'grey' ) 







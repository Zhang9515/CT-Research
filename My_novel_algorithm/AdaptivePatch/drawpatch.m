clear all
load ..\..\..\Data\Adaptive_patchsize_selection\trial2D_prior_324
pic = trial2D_prior ;
figure
set (gcf,'Position',[0,0,512,512]);   
for i = 1 :3
    for j =1 : 3
         imshow(pic(317-1+i-40:317-1+i+40,160-1+j-40:160-1+j+40),[0 0.5],'border', 'tight','initialmagnification','fit');       
         add = ['..\..\..\Data\Adaptive_patchsize_selection\', num2str(i+(j-1)*3)  ];
         eval(['print -djpeg -r600 ', add]);
    end
end

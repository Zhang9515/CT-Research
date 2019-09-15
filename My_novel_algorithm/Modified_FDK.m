% 2017/04/13 ZXZ
%%  load image

a = open ( 'D:\TestCpp\CT\Data\FDK\Modified_Para\saddle\For paper\H5Dis103middleP0.2.fig' ) ;
h = get (gca , 'Children') ;
data = get ( h , 'Cdata') ;
% data = flipud ( data ) ;
PicSize = 513 ;

Size = [ 60 , 60 , 60 ] ;     % actual range
Resolution2 = max ( Size ) / PicSize ; 
Center_t = max ( Size ) / 2 ;  Center_s = max ( Size ) / 2 ;   Center_z = max ( Size ) / 2 ;          % define the center 

%% modification   circular trajectary
% modifying function :  Err = Intensity * ( Z / D )^2 * ( K * r^2 + C ) ;
% Distance = 103 ;          % distance between source and center point
% MFDK = zeros ( PicSize , PicSize ) ;
% K = 20 ; 
% C = 3.5e3 ; 
% Tindex = 257 ;
% zm = max ( Size ) / 2 ; rm = max ( Size ) / 2 ;
% for Zindex = 1 : PicSize
%         for Sindex = 1 : PicSize
%                 z = ( Zindex - 0.5 ) * Resolution2 - Center_z ; s = ( Sindex - 0.5 ) * Resolution2 - Center_s ; 
%                 t = ( Tindex - 0.5 ) * Resolution2 - Center_t ;
%                 r = sqrt ( s^2 + t^2 ) ; 
%                 MFDK ( Zindex , Sindex ) =  data ( Zindex , Sindex ) / ( 1 - ( z / ( zm * Distance ) )^2 * ( K * ( r / rm )^2 + C ) ) ;
%         end
% end
% figure , imshow ( MFDK(102:412,:) , [ 0.15,0.35 ] ) ;
% figure , plot (   1 : 512 , squeeze ( data0 ( : , 256  ) ) ,1 : 512 , squeeze ( MFDK ( : , 256  ) ) , 1 : 512 , squeeze ( data ( : , 256  ) ) ) ; axis ( [ 102 412 0.1,0.35 ] ) ;
% set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;     
% set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;
% figure , plot (   1 : 512 , squeeze ( data0 ( : , 257  ) ) ,1 : 512 , squeeze ( data ( : , 257  ) ) ) ; axis ( [ 102 412 0.15,0.35 ] ) ;
% set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;     
% set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;
% figure , plot (   1 : 512, squeeze ( data ( : , 256 ) ) ,1 : 512 , squeeze ( MFDK ( : , 256 ) ) , 1:512 , squeeze ( SFBPimage ( 256 , 256 , : ) ) ) ;
% set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;     
% 296 low-density
% set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;axis ( [ 102 412 0.4 1.1 ] ) ;
%% modification saddle-line trajectary
% modifying function : Err  = Intensity * ( C6 * ( z -  ( C1  + C2 * R ^ 2 ) * H ^ 2 * sin ( 2B ) ) ^ 2 + H ^ 2 * ( C3 + C5 * sin ( 2B ) ) ) / ( D^2 )
% 
%  
Distance = 103.9 ;          % distance between source and center point
H  = 5 ;                     % height of saddle line
Beta = 0;            % angle of saddle line 
MFDK = zeros ( PicSize , PicSize ) ;

C3 = 0.5 ; 
C6 = 0.4 ;
Tindex = 257 ;
for Zindex = 1 : PicSize
        for Sindex = 1 : PicSize
                z = ( Zindex - 0.5 ) * Resolution2 - Center_z ; s = ( Sindex - 0.5 ) * Resolution2 - Center_s ; 
                t = ( Tindex - 0.5 ) * Resolution2 - Center_t ;
                r = sqrt ( s^2 + t^2 ) ;
                MFDK ( Zindex , Sindex ) =  data ( Zindex , Sindex ) / ( 1 - ( C6 * ( z ) ^ 2 + H^2 * ( C3 ) ) / ( Distance^2 ) ) ;
        end
end

figure , imshow ( MFDK ( 102 : 412 , : ) , [ 1,1.05 ] ) ;
figure , plot ( 1 : 513 , squeeze ( data0 ( : , 296 ) ) , 1 : 513 , MFDK ( : , 296 ) ) ;
set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;
set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;axis ( [ 102 412 0.98 1.05 ] ) ;
figure , plot ( 1 : 513 , squeeze ( data0 ( : , 296 ) ) , 1 : 513, squeeze ( data ( : , 296 ) ) ) ;
set ( gca , 'XTick' , [ 0 , 56 , 156 , 256 , 356 , 456 , 512 ] ) ;
set ( gca , 'XTickLabel' , { '256' , '200' , '100' , '0' , '-100' , '-200' , '-256' } ) ;axis ( [ 102 412 0.98 1.05 ] ) ;
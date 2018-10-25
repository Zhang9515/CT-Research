Display = zeros ( height , width ) ;        % store reconstruction result
for j = 1 : height                                     % rebinning
    for i = 1 : width
%             i= 90 ; j =101;
            X = ( i - 0.5 ) * Resolution2 - Center_x ; Y = ( j - 0.5 ) * Resolution2 - Center_y ; 
            Xindex = round ( ( X - XSample1 ) / ResolutionPI ) + 1 ; Xindex1 = floor ( ( X - XSample1 ) / ResolutionPI ) + 1 ; Xindex2 = Xindex1 + 1 ;                 % location 0.5 represents pixel density
            x1 = X1 + ( Xindex1 - 0.5 ) * ResolutionPI ; x2 = X1 + ( Xindex2 -0.5 ) * ResolutionPI ;
            Yindex = round ( ( Y - YSample1 ) / ResolutionPI ) + 1 ; Yindex1 = floor ( ( Y - YSample1 ) / ResolutionPI ) + 1 ; Yindex2 = Yindex1 + 1 ; 
            y1 = Y1 + ( Yindex1 - 0.5 )* ResolutionPI ; y2 = Y1 + ( Yindex2 - 0.5 )* ResolutionPI ;
            if ( Xindex1 >= 1 && Yindex1 >= 1 &&  Xindex2 <= PIWidth && Yindex2 <= PIHeight ) 
                 Display ( j , i ) = ( ( y2 - Y ) * ( x2 - X ) * PISpaceCt ( Yindex1 , Xindex1 ) + ( y2 - Y ) * ( X - x1 ) * PISpaceCt ( Yindex1 , Xindex2 ) ...
                     + ( Y - y1 ) * ( x2 - X ) * PISpaceCt ( Yindex2 , Xindex1 ) + ( Y - y1 ) * ( X - x1 ) * PISpaceCt ( Yindex2 , Xindex2 ) ) / ResolutionPI^2 ;
            end
    end    
end
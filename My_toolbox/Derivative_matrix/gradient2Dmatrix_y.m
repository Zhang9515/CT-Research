function Gradienty = gradient2Dmatrix_y( height , width )
    % sobel operator, using mirror extending
%     Gradienty = zeros( height * width ) ;
    Sobel = [ -1 -2 -1 1 2 1 ] ;
    gyRow=[]; gyCol=[]; gyW=[]; 
    for x = 1 : width
        for y = 1 : height
             gyRow(end+(1:6)) = (y-1) * height + x ;
            if ( x == 1)
                xindex = [0,0,1] ;
            elseif ( x == width)
                xindex = [-1,0,0];
            else 
                xindex = [-1,0,1];
            end
            if ( y == 1)
                yindex = [0,1] ;
            elseif ( y == height)
                yindex = [-1,0];
            else
                yindex = [-1,1];
            end
            [xindex,yindex] = meshgrid(xindex,yindex);
            xindex = xindex' ; yindex = yindex' ; 
             % result is a matrix in a x-y cordinate but here i need extract
            % them in linear way in which x-first then-y, because matlab
            % extract data by column (y), so i transpose the matrix.
            for n =1 : 6
                gyCol(end+1) = ((y + yindex(n)) -1) * width + (x + xindex(n)) ;
                gyW(end+1) =  Sobel(n) ;
            end 
        end
    end
    Gradienty = sparse ( gyRow , gyCol , gyW , height * width ,  height * width ) ;
end
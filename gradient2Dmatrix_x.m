function Gradientx = gradient2Dmatrix_x( height , width )
    % sobel operator, using mirror extending
    Sobel = [ -1 1 -2 2 -1 1 ] ;
    gxRow=[]; gxCol=[]; gxW=[]; 
    for x = 1 : width
        for y = 1 : height
            gxRow(end+(1:6)) = (y-1) * height + x ;
            if ( x == 1)
                xindex = [0,1] ;
            elseif ( x == width)
                xindex = [-1,0];
            else
                xindex = [-1,1];
            end
            if ( y == 1)
                yindex = [0,0,1] ;
            elseif ( y == height)
                yindex = [-1,0,0];
            else
                yindex = [-1,0,1];
            end
            [xindex,yindex] = meshgrid(xindex,yindex);
            % result is a matrix in a x-y cordinate but here i need extract
            % them in linear way in which x-first then-y, because matlab
            % extract data by column (y), so i transpose the matrix.
            xindex = xindex' ; yindex = yindex' ; 
            for n =1 : 6
                gxCol(end+1) = ((y + yindex(n)) -1) * width + (x + xindex(n)) ;
                gxW(end+1) =  Sobel(n) ;
            end 
        end
    end
    Gradientx = sparse ( gxRow , gxCol , gxW , height * width ,  height * width ) ;
end
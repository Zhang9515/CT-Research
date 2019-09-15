function divergence_matrix = divergenceMatrix2D( height ,width ) 
    Gradientx = gradient2Dmatrix_x( height , width ) ;
    Gradienty = gradient2Dmatrix_y( height , width ) ;
    divergence_matrix = Gradientx'*Gradientx + Gradienty'*Gradienty ;
end
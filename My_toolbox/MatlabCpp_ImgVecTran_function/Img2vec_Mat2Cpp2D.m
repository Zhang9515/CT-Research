function picvector = Img2vec_Mat2Cpp2D( Img )
    picvector = single( reshape ( fliplr(Img'), numel(Img) , 1) ) ;
end
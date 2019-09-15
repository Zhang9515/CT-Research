function X = col_normalization ( inputmatrix )
    [Rnum , Cnum] = size(inputmatrix) ; 
    X = zeros(Rnum , 1) ;
    for i = 1 : Cnum
        if ( abs (norm(inputmatrix(:,i)) ) > eps )
            X( : , end + 1 ) = inputmatrix(:,i) / ( norm(inputmatrix(:,i))) ;
        end
    end
    X( : , 1) =[] ;
end
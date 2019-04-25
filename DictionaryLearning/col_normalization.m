function X = col_normalization ( inputmatrix )
    [Rnum , Cnum] = size(inputmatrix) ; 
    X = zeros(Rnum , 1) ;
    eps = 10^-10 ;   % to avoid zero in the denominator
    for i = 1 : Cnum
        if ( abs (norm(inputmatrix(:,i)) ) > eps )
            X( : , end + 1 ) = inputmatrix(:,i) / ( norm(inputmatrix(:,i))) ;
        end
    end
    X( : , 1) =[] ;
end
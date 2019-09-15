function  Result = isotropicISTA ( threshold , normalizer)
% X : from left to right; Y : from down to up; Z : from near to far
% input: image: N(=X-Y-Z)-3, ps: length of patch; Result:  N( =X*Y*Z )-3 
    EPS = 1e-15 ;     
    sizeiter = size(threshold ,1 ) ; 
    Result = zeros( sizeiter , 3 ) ;
    for num = 1 : sizeiter
        vector = threshold( num , : ) ; 
        Result(num , : ) = vector * max( norm( vector , 2 ) - normalizer , 0 ) / (norm( vector , 2 ) + EPS) ;
    end

end
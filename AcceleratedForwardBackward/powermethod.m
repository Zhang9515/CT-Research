function  [l , s ] = powermethod( A , x0 , eps )
% power method to get the max eigenvalue of the matrix
    if nargin==2
        eps = 1.0e-6 ;    
    end
    v = x0;    % arbitary vector
    M = 5000;  % max bound of iterate time 
    
    m = 0;
    l = 0;
    for ( k = 1 : M )
        y = A * v ;
        m = max( y ) ; % m is the largest component of y
        v = y / m ;
        if ( abs( m - l ) < eps )
            l = m ; % stop  and l is max eigen value
            s = k ; % s is iterate time
            return ;
        else
            if ( k == M )
                disp ( 'Reach the max iterate time !' ) ;
                Err = abs( m - l ) ; 
                disp ( Err ) ;
                l = m ;
                s = M ;
            else
                l = m ;
            end
        end
    end

return;
function [delta_x , delta_y , delta_wx , delta_wy , wtv ] = WTV ( U , height , width ) 
% weighted TV, U is vector, diretion : 1 is x ; 2 is y   
        
    H = []; W = []; V = [];
    H( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width +1;
    W( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width +1;
    V( (end + 1) : (end + height) ) = -1;

    H( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width +1;
    W( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width +2;
    V( (end + 1) : (end + height) ) = 1;

    for i = 1 : height
        H( (end + 1) : (end + (width-2)) ) =  ( i -1) * width + (2 : (width-1)) ;
        W( (end + 1) : (end + (width-2)) ) = ( i -1) * width + (1 : (width-2)) ; 
        V( (end + 1) : (end + (width-2)) ) = -0.5 ;

        H( (end + 1) : (end + (width-2)) ) = ( i -1) * width + (2 : (width-1)) ;
        W( (end + 1) : (end + (width-2)) ) = ( i -1) * width + (3 : width) ;
        V( (end + 1) : (end + (width-2)) ) = 0.5 ;
    end

    H( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width + width ;
    W( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width + width - 1 ;
    V( (end + 1) : (end + height) ) = -1 ;

    H( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width + width ;
    W( (end + 1) : (end + height) ) = ( (1 : height) - 1) * width + width ;
    V( (end + 1) : (end + height) ) = 1 ;

    delta_x = sparse ( H , W , V , height * width , height * width ) ;  % build the sparse matrix as gradient matrix
       
    H = [] ; W = [] ; V = [] ;

    H( (end + 1) : (end + width) ) = 1 : width ;
    W( (end + 1) : (end + width) ) = 1 : width ;
    V( (end + 1) : (end + width) ) = -1;

    H( (end + 1) : (end + width) ) = 1 : width ;
    W( (end + 1) : (end + width) ) = (1 : width) + width ;
    V( (end + 1) : (end + width) ) = 1;

    for i = 1 : width
        H( (end + 1) : (end + (height-2)) ) =  ( (2 : ( height - 1 )) - 1 ) * width + i ;
        W( (end + 1) : (end + (height-2)) ) = ( (1 : ( height - 2 )) - 1 ) * width + i ; 
        V( (end + 1) : (end + (height-2)) ) = -0.5 ;

        H( (end + 1) : (end + (height-2)) ) = ( (2 : ( height - 1 )) - 1 ) * width + i ;
        W( (end + 1) : (end + (height-2)) ) = ( (3 : height) - 1 ) * width + i ;
        V( (end + 1) : (end + (height-2)) ) = 0.5 ;
    end

    H( (end + 1) : (end + width) ) = ( height - 1) * width + (1 : width) ;
    W( (end + 1) : (end + width) ) = ( height - 2) * width + (1 : width) ;
    V( (end + 1) : (end + width) ) = -1 ;

    H( (end + 1) : (end + width) ) = ( height - 1) * width + (1 : width) ;
    W( (end + 1) : (end + width) ) = ( height - 1) * width + (1 : width) ;
    V( (end + 1) : (end + width) ) = 1 ;
        
    delta_y = sparse ( H , W , V , height * width , height * width ) ;
    
    Display_x = delta_x * U ; 
    Display_y = delta_y * U ; 
    
    W_x = WTV_weight ( Display_x ) ;
    W_y = WTV_weight ( Display_y ) ;
    
    delta_wx = repmat( W_x , 1 , height * width ) .* delta_x ; 
    delta_wy = repmat( W_y , 1 , height * width ) .* delta_y ;     
    
    wtv = norm( W_x .* Display_x , 1 ) + norm( W_x .* Display_x , 1 ) ;

return;

function wtv_weight = WTV_weight ( tv , miu )

if nargin == 1
    miu = 15 ;
end
wtv_weight = 1 ./ ( ( miu * log( 2 ) ) .* ( 1 + exp( abs (tv) ./ miu ) ) ) ;

return;
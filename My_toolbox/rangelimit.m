function output = rangelimit ( input )
if ( input > 1e16 )
    output = 1e16 ;
elseif ( input < -1e16 )
    output = -1e16 ;
else
    output = input ;
end
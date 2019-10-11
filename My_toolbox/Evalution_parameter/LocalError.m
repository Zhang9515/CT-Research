function local_error = LocalError( current , previous )
    local_error = norm( current - previous , 2) / norm(previous,2) ;
end
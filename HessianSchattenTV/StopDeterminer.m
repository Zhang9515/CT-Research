function result = StopDeterminer( VarCur , VarPre , threshold )
% use local relative error
    datalength = numel( VarCur) ; 
    VarCur = reshape ( VarCur , datalength , 1 ) ;
    VarPre = reshape ( VarPre , datalength , 1 ) ;
    
    LocRelErr = norm( VarCur - VarPre , 2 ) / ( norm( VarPre , 2 ) + 1e-10 ) ; 
    if ( LocRelErr <= threshold )
            result = true ;
    else
            result = false ;
    end
    disp ( ['local relative error: ', num2str( LocRelErr ) ] ) 
end
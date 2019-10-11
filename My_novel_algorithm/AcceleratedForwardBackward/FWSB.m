function U = FWSB ( delta_wx , delta_wy , V , alpha , lamda , step, MinLim , MaxLim  )
% fast weighted spliting bregman method  
    
    U0 = V ;
    Uprevious = U0 ;
    ex0 = 0; ey0 = 0;
    ex_previous = ex0 ; ey_previous = ey0 ;
   
    for jLoop = 1 : 10
            
            Ux = delta_wx * Uprevious ; 
            Uy = delta_wy * Uprevious ;
%             Xx = Ux ; Xy = Uy;     % initial inner variable
            Xx = zeros(size(Ux)) ; Xy = zeros(size(Uy)) ;
            Zx = Ux + ex_previous ; Zy = Uy + ey_previous ; 
            ex = Cut( Zx , alpha / (lamda+eps) ) ; ey = Cut( Zy , alpha / (lamda+eps) ) ; 
            
            Xprevious = Uprevious ; 
            for mLoop = 1 : 100                
                X = V - step * lamda * ( delta_wx' * ( Xx + 2 * ex - Zx ) + delta_wy' * ( Xy + 2 * ey - Zy ) ) ;
                X( X>MaxLim ) = MaxLim ; X( X<MinLim ) = MinLim ; 
                local_error = LocalError( X,Xprevious ) ;      
                Xprevious = X ; 
                Xx = delta_wx * Xprevious ; Xy = delta_wy * Xprevious ;  
                residual = norm(X - V + step * lamda * ( delta_wx' * ( Xx + 2 * ex - Zx ) + delta_wy' * ( Xy + 2 * ey - Zy ) ) ) ; 
                disp( [ '     mLoop: ' , num2str( mLoop ) , '   local_e : ' , num2str( local_error ), ' residual: ' , num2str(residual)] ) ;
            end      %mLoop
            U = X ; 
            loss = 0.5 / step * norm( U-V ) + alpha * (norm(delta_wx * U,1) + norm(delta_wy * U,1)) ;
            disp( [ '   jLoop: ' , num2str( jLoop ) , '   local_e : ' , num2str( LocalError( U,Uprevious ) ) , ' loss: ', num2str( loss ) ] ) ;
            Uprevious = U ;
    end % jLoop

end % function
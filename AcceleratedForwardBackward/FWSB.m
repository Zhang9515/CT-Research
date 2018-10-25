function X = FWSB ( delta_wx , delta_wy , V , lamda , theta , beta )
% fast weighted spliting bregman method  
    
    U0 = V ;
    Uprevious = U0 ;
    ex0 = 0; ey0 = 0;
    ex_previous = ex0 ; ey_previous = ey0 ;
   
    jLoop = 1 ;
    while(1)
        disp( [ '   jLoop:' , num2str( jLoop ) ] ) ;
        Ux = delta_wx * Uprevious ; 
        Uy = delta_wy * Uprevious ;
        Xx = Ux ; Xy = Uy;
        Zx = Ux + ex_previous ; Zy = Uy + ey_previous ; 
        ex = Cut( Zx , lamda / theta ) ; ey = Cut( Zy , lamda / theta ) ; 
        
        mLoop = 1 ;
        Xprevious = Uprevious ; 
        while(1)
            disp( [ '     mLoop:' , num2str( mLoop ) ] ) ;
            X = V - beta * theta * ( delta_wx' * ( Xx + 2 * ex - Zx ) + delta_wy' * ( Xy + 2 * ey - Zy ) ) ;
            if ( norm ( Xprevious - X , 2) <= 10e-6 * norm( Xprevious , 2 ) )
                break;
            end
            Xprevious = X ; 
            Xx = delta_wx * Xprevious ; Xy = delta_wy * Xprevious ;  
            mLoop = mLoop + 1 ;
        end
        
        if ( norm ( Uprevious - X , 2) <= 10e-6 * norm( Uprevious , 2 ) )
                break;
        end
        Uprevious = X ; 
        jLoop = jLoop + 1 ;
    end

return;

function output = Cut ( input, threshold )
    output = zeros(length(input),1) ;
    for i = 1 : length(input)
        if  abs(input( i )) < threshold
              output(i) = input( i ) ;
        elseif  input( i ) > threshold
              output(i) = threshold ;  
        elseif input( i ) < (-threshold)
              output(i) = - threshold ;  
        end
    end

return;
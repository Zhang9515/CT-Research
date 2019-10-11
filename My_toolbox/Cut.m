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
end
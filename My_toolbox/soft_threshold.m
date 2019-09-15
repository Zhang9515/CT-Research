function [ x ] = soft_threshold( b,T )
    x = sign(b).*max(abs(b) - T,0);
end

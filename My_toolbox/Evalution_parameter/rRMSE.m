function value = rRMSE(input,ref)
    value = sqrt(mean2(input - ref)) / (sqrt(mean2(ref))+eps) ; 
end
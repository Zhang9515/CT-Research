function value = rRMSE(input,ref)
    value = sqrt(mean2((input - ref).^2)) / (sqrt(mean2(ref.^2))+eps) ; 
end
function rmse = RMSE( result , ref)
    rmse = sqrt( mean2((result - ref).^2) );
return
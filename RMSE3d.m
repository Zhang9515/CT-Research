function rmse = RMSE3d(result,ref)
rmse = sqrt(mean(mean(mean((result-ref).^2)))) / numel(result);
return
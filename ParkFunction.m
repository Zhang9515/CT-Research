% Park_Window for shortscan
% 2018/03/23 by ZXZ
% all angle should be radian 
% MaxGama is the Fan beam angle not half  Fan beam angle
% output is the weight of the specific projection line
function Weight = ParkFunction ( beta , betaStartAngle , gama , MaxGama )

if ( beta >= betaStartAngle && beta < betaStartAngle + MaxGama - 2 * gama )
    Weight = sin ( pi / 4 * beta / ( MaxGama / 2 - gama ) )^2 ;
elseif ( beta >= betaStartAngle + MaxGama - 2 * gama && beta <= betaStartAngle + pi - 2 * gama )
    Weight =  1 ;
elseif ( beta > betaStartAngle + pi - 2 * gama && beta <= betaStartAngle + pi + MaxGama )
    Weight = sin ( pi / 4 * ( pi + MaxGama - beta ) / ( MaxGama / 2 + gama ) )^2 ;
else 
    Weight = 0 ;
end
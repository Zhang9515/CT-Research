function Fai = ObjectiveFunction( ProjectionData, AU , U, OrderofHS, Tao)
    Hsnorm
    Fai = 0.5 * (ProjectionData - AU)' * (ProjectionData - AU) + Tao * Hsnorm ;

end
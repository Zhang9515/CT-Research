function solution = Cholesky( L , d)
        y(1)=d(1);
        n = length(d);       
        for i=2:n   
            for j=1:i-1
                d(i)=d(i)-L(i,j)*y(j);
            end
            y(i)=d(i);
        end
        U = L';
        x(n)=y(n)/U(n,n);
        for i=(n-1):-1:1
            for j=n:-1:i+1
                y(i)=y(i)-U(i,j)*x(j);
            end
            x(i)=y(i)/U(i,i);    
        end
        solution = x';
end
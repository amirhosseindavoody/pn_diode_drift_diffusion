function J = LU_decomposition (A,f)
    [n,m] = size(f); % n = number of total space meshes
    %--------------------Setting Values of The Upper and Lower Trangular
    %--------------------Matrix According to A Matrix Values
    for i = 1:n
        if(i == 1)
            L(i,i) = 1;
            U(i,i) = A(i,i);
            U(i,i+1) = A(i,i+1);
        elseif (i < n)
            L(i,i) = 1;
            L(i,i-1) = A(i,i-1)/U(i-1,i-1);
            U(i,i) = A(i,i)-L(i,i-1)*A(i-1,i);
            U(i,i+1) = A(i,i+1);
        else
            L(i,i) = 1;
            L(i,i-1) = A(i,i-1)/U(i-1,i-1);
            U(i,i) = A(i,i)-L(i,i-1)*A(i-1,i);
        end;
    end;
    
    %--------------------Solving the equation L*y = f
    for i = 1:n
        m = f(i,1);
        for j = 1:i
            if (j<i)
                m = m-L(i,j)*y(j,1);
            else
                y(j,1) = m/L(i,j);
            end;
        end;
    end;
    
    %-----------------------------Solving the equation U*V = y
    for i = 1:n
        m = y(n-i+1,1);
        for j = 1:i
            if (j<i)
                m = m-U(n-i+1,n-j+1)*J(n-j+1,1);
            else
                J(n-j+1,1) = m/U(n-i+1,n-j+1);
            end;
        end;
    end;
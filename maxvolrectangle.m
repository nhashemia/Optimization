function [L, U, phi, Hessian] = maxvolrectangle(A,b)

[m,n] = size(A);
r_A = rank(A);

if r_A ~= n
    error('Inapproperiate A');
else
    no_c = m; % Number of constraints
end

% Number of constraints are greater than the dimension of the space m>n.
A_r = A';
[~,basiccol] = rref(A_r);
A_r = A_r(:,basiccol)';
b_r = b(basiccol,1);

[r_A,c_A] = size(A);
A_plus = zeros(r_A,c_A);
A_minus = zeros(r_A,c_A);

for i = 1:r_A
    for j = 1:c_A
        if A(i,j) > 0
            A_plus (i,j) = A(i,j);
            A_minus (i,j) = 0;
        else
            A_plus (i,j) = 0;
            A_minus (i,j) = -A(i,j);
        end        
    end
end

% A_r : independent rows of A.
% rank(A) = n.
% As the vertices of the polytope are the possible optimal solutions for
% the vertices of the rectangle:

% V_R = inv(A_r)*b_r; 
% L(:,1) = V_R - V_R*0.5;
% U(:,1) = V_R ;
L(1:n,1) = 0.1*ones(1:n,1);
U(1:n,1) = 0.2*ones(1:n,1);

miu = 1.5;
epsilon = 0.01;
delta = zeros(1:2*n,1);
Hessian = zeros(1:2*n,1:2*n);
counter = 1;
k(1) = 0.5;
phi(1) = 0;
g = zeros(1:r_A,1);

while m/k(counter) >= epsilon
    for i = 1:r_A
        g(i) = A_plus(i,:)*U(:,counter) - A_minus(i,:)*L(:,counter) - b(i);  
        if g(i) > 0
            phi (counter) = 100000; % the point is out of the domain.
            i = r_A + 1; 
        else
            phi (counter+1) = phi (counter) - log10(-g(i));
        end
    end
    
    if phi (counter) == 100000
        counter = counter + 1;
        X_temp = [U(:,counter-1);L(:,counter-1)] - 0.1.*inv(Hessian)*delta(:,counter-1);
        U(:,counter) = X_temp(1:n);
        L(:,counter) = X_temp(n+1:2*n);
        k(counter) = k(counter-1)*miu;   
    else
        for j = 1:n
            delta(j,counter) = -k(counter)/(U(j,counter) - L(j,counter)) - (1./g)*A_plus(:,j);
            delta(j+n,counter) = k(counter)/(U(j,counter) - L(j,counter)) + (1./g)*A_minus(:,j);
            for t = 1:n
                if j == t
                    Hessian(j,j) = k(counter)/((U(j,counter) - L(j,counter))^2) + (1./(g.*g))*(A_plus(:,j).^2);
                    Hessian(j+n,j+n) = k(counter)/((U(j,counter) - L(j,counter))^2) + (1./(g.*g))*(A_minus(:,j).^2);
                    
                    Hessian(j,t+n) = -k(counter)/((U(j,counter) - L(j,counter))^2) - (1./(g.*g))*(A_plus(:,j).*A_minus(:,t));
                    Hessian(j+n,t) = -k(counter)/((U(j,counter) - L(j,counter))^2) + (1./(g.*g))*(A_minus(:,j).*A_minus(:,t));
                else
                    Hessian(j,t) = (1./(g.*g))*(A_plus(:,j).*A_plus(:,t));
                    Hessian(j,t+n) = -(1./(g.*g))*(A_plus(:,j).*A_minus(:,t));
                    Hessian(j+n,t) = (1./(g.*g))*(A_minus(:,j).*A_minus(:,t));
                    Hessian(j+n,t+n) = -(1./(g.*g))*(A_minus(:,j).*A_plus(:,t));
                end
            end
        end    
        counter = counter + 1;
        X_temp = [U(:,counter-1);L(:,counter-1)] - 0.1.*inv(Hessian)*delta(:,counter-1);
        U(:,counter) = X_temp(1:n);
        L(:,counter) = X_temp(n+1:2*n);
        k(counter) = k(counter-1)*miu;   
    end
end
end


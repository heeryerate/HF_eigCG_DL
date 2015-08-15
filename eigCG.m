function [V, T, sol] = eigCG(A, b, P, x, m, nev, TOL, MAXITER)
%% Initialization

assert(size(A, 1) > m, 'm should be smaller than n!')
assert(m > 2*nev, 'm should be greater than 2*nev!')

r0 = b - A * x;
r = r0;
z = inv(P) * r;
p = z;
rho_pre = r' * r;

alpha = 0;
beta = 0;
flag = 1;

i = 1;
j = 1;
%% CG iteration
while(1)
    
    %     Step 10.1 Can be deleted
    %     if i == m
    %         t_pre = t;
    %     end
    
    t = A *p;
    alpha1 = alpha;
    alpha = r' * z / (p' * t);
    
    % Initialize V be set of Lanczos vectors
    V(:, i) = z / sqrt(rho_pre);           
    
    x = x + alpha * p;
    r = r - alpha * t;
    z = inv(P) * r;
    rho = r' * z;
    beta1 = beta;
    rho_pre1 = rho_pre;
    beta = rho / rho_pre;
    rho_pre = rho;
    p = z + beta * p;
    
    % Step 11 eigCG additions
    % Generate Tridiagonalize matrix T from CG procedure
    if i == 1
        T(i, i) = 1 / alpha;
    elseif (i > 1)
        T(i, i) = 1 / alpha + beta1 / alpha1;
        T(i-1, i) = -sqrt(beta1) / alpha1;
        T(i, i-1) = T(i-1, i);
    end
    
    %     if i < m && flag == 1
    %          fprintf('norm(V'' * A * V - T) = %f: \t Tridiagonalize matrix A to T in %d dimensional space:\n', norm(V' * A * V - T), i)
    %     end
    
    % flag: indicator of RESTART, to make sure RESTART happens exactly once
    if i == m && flag == 1    
        Tm = T(1:m, 1:m);
        [Y, M] = eigs(Tm, nev, 'sm');
        [Y1, M1] = eigs(T(1:(m-1), 1:(m-1)), nev, 'sm');
        Y1 = [Y1; zeros(1, nev)];
        rankbuf = rank([Y, Y1]);
        if rankbuf ~= 2* nev
            warning('warning: [Y, Y1] are not rank 2*nev')
        end
        Q = orth([Y, Y1]);
        H = Q' * Tm * Q;
        [Z, M] = eig(H);
        
        % RESTART
        V = V * Q * Z;
        i = rankbuf;
        T(1:i, 1:i) = M;
        %         w =  t - beta1 * t_pre;
        %         wbuff = w' * V / rho_pre;
        %         wbuff1 = wbuff(1:i);
        %         T(i+1, 1:i) = wbuff1;
        %         T(1:i, i+1) = wbuff1';
        flag = 0;
        fprintf('RESTART happens at CG step j = %d.\n', j)
    end
    
    % Terminating condition
    if (norm(r)/norm(r0) < TOL || j > MAXITER) && flag == 0
        break;
    end
    
    i = i + 1;
    j = j + 1;
end
%% Output
sol = x';
% if j < m
%     warning('warning: no restart!')
% end
fprintf('Relative error is %f at CG step j = %d.', norm(r)/norm(r0), j)
end

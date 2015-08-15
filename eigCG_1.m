function [sol, vecj, err, v, valbk] = eigCG(Afunc, b, x0, maxiters, miniters, Mdiag)
%function [V, T, sol] = eigCG(A, b, P, x, m, nev, TOL, MAXITER)
%% Initialization

% assert(size(A, 1) > m, 'm should be smaller than n!')
% assert(m > 2*nev, 'm should be greater than 2*nev!')

MAXITER = maxiters;
x = x0;
TOL = 5e-4;
eps = 1e-8;
inext = 5;
imult = 1.3;
vecj = [];
sol = {};

gapratio = 0.1;
mingap = 10;

maxtestgap = max(ceil(maxiters * gapratio), mingap) + 1;

vals = zeros(maxtestgap,1);


Ax = Afunc(x);
nev = 8;
m = min(ceil(size(Ax, 1) * 0.05), 40);


r0 = b - Ax;
r = r0;
z = r ./ Mdiag;
p = z;
rho_pre = r' * r;

alpha = 0;
beta = 0;
flag = 1;

i = 1;
j = 1;
normb = norm(b);
val = 0.5*double((-b-r)'*x);
valbk = 0;
v = 0;
%% CG iteration
while(1)
    
    %     Step 10.1 Can be deleted
    %     if i == m
    %         t_pre = t;
    %     end
    
    t = Afunc(p);
    alpha1 = alpha;
    alpha = r' * z / (p' * t);
    
    % Initialize V be set of Lanczos vectors
    V(:, i) = z / sqrt(rho_pre);
    
    x = x + alpha * p;
    
    r = r - alpha * t;
    z = r ./ Mdiag;
    rho = r' * z;
    beta1 = beta;
    rho_pre1 = rho_pre;
    beta = rho / rho_pre;
    rho_pre = rho;
    p = z + beta * p;
    
    val = 0.5*double((-b-r)'*x);
    vals( mod(j-1, maxtestgap)+1 ) = val;
    testgap = max(ceil( j * gapratio ), mingap);
    prevval = vals( mod(j-testgap-1, maxtestgap)+1 ); %testgap steps ago
    
    if j == ceil(inext)
        vecj(end+1) = j;
        sol{end+1} = x;
        inext = inext*imult;
    end
    
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
        %i = rankbuf;
        i = size(Q, 2);
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
    %    err = norm(r)/norm(r0); && prevval < 0
    
    err = (val - prevval)/val;
    termi = (j > testgap && prevval < 0 && err < TOL*testgap && j >= miniters);
    
    if flag == 0
        T = T(1:nev, 1:nev);
        eigT = sort(diag(T)')
        if eigT(1) < - eps * normb
            warning('warning: Negative curvature captured!')
            v = V(: ,1);
            if v' * b >= 0
                v = 0;
                valbk = 0;
                break;
            end
            valbk = -0.5*double((-b-r)'*v);
            warning('This is a negative curvature direction')
            break;
        end
        flag = 2;
    end
    
    if normb <= 1
        if termi && (flag == 2)
            break;
        end
    else
        if termi
            break;
        end
    end
    
    
    
    i = i + 1;
    j = j + 1;
end
%% Output
% sol = x';
% if j < m
%     warning('warning: no restart!')
% end

fprintf('**************norm(grad) = %f\n', normb);
fprintf('**************norm(Ax-b) = %f\n', norm(Afunc(x)-b));
if j ~= ceil(inext)
    vecj(end+1) = j;
    sol{end+1} = x;
end
fprintf('Relative error is %f at CG step j = %d.\n', err, j)
end

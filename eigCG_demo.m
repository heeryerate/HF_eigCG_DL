clear
clc
rng(1);

n = 500;
n = size(A, 2)
m = 130;
nev = 8;

%% Generate dense matrix
% A = rand(n, n);
% A = A' + A;
% mineigA = min(neig(A, nev));
% if mineigA < 0
%     A = A -  (mineigA - 1) * eye(n);
% end
% assert(min(neig(A, nev)) > 0, 'error: A should be positive definite')

%% Generate sparse matrix
% density = 0.1;
% rc = 5 * rand(n, 1) + 0 * ones(n, 1);
% A = sprandsym(n, density, rc);

%% eigCG
b = A * ones(n, 1);
P = eye(n);
x = rand(n,1)*10;
TOL = 1e-06;
MAXITER = 1e+3;

[V, T, sol] = eigCG(A, b, P, x, m, nev, TOL, MAXITER);

%% Output
V = V(:, 1:nev);
T = T(1:nev, 1:nev);
% V' * A * V - T 
fprintf('\nnev lowest eigenvalues of A:\n')
eigA = sort(rc);
%eigA = neig(A, nev)';
fprintf('\t%f', eigA(1: nev))
eigT = sort(diag(T)');
fprintf('\nnev lowest eigenvalues by eigCG:\n')
fprintf('\t%f', eigT(1: nev))
fprintf('\n')
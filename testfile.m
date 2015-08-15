% test file
clear
clc
rng('default');
rng(1111)

n = 1000;
m = 40;
nev = 8;

density = 0.1;
rc = sign(rand(n, 1) - 0.5) .* rand(n, 1);
%rc = 5 * rand(n, 1) + 0 * ones(n, 1);
global A;
A = sprandsym(n, density, rc);

%% eigCG
%b = A * ones(n, 1);
b = rand(n,1)/100;
P = eye(n);
x = rand(n,1)*10;

tic
[xs, is, err] = conjgrad_2(@foo, b, x, 1000, 1, 1);
toc
%is
norm(A * xs{end} - b)
% a = norm(xs{end});

tic
min(rc)
[xs, is, err, v, valbk] = eigCG_1(@foo, b, x, 100, 1, 1);
%[xs, is] = conjgrad_2(@foo, b, x, 1000, 1, 1);
toc
%is
norm(A * xs{end} - b)
% b = norm(xs{end});
% a - b
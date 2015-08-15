function afun = foo(x)
% A = diag(1:10)./10;
global A;
assert(size(A, 2) == size(x, 1), 'Error: dimension not match')
afun = A * x;
end
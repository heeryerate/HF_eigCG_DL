function eigA = neig(A, nev)
% Get the lowest nev eigenvalues
eigA = eig(A);
eigA = eigA(1:nev);
end
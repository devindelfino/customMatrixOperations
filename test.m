% Filename: test.m
% Description: Testing environment to validate custom operation implementations against the corresponding built-in MATLAB functions.

% ------------------------------------------------

Am = randi([-15 15], 7);
Bm = randi([-15 15], 7);
Ac = Matrix(Am);
Bc = Matrix(Bm);

assert(all((Ac+Bc) == (Am+Bm)), 'Matrix Addition FAILED.');
display('Testing Matrix Addition successful.')

% ------------------------------------------------

Am = randi([-15 15], 3);
Bm = randi([-15 15], 3);
Ac = Matrix(Am);
Bc = Matrix(Bm);

assert(all((Ac-Bc) == (Am-Bm)), 'Matrix Subtraction FAILED.');
display('Matrix Subtraction successful.')

% ------------------------------------------------

Am1 = randi([-15 15], 7);
Bm1 = randi([-15 15], 7);
Ac1 = Matrix(Am1);
Bc1 = Matrix(Bm1);

Am2 = randi([-15 15], 3);
Bm2 = randi([-15 15], 3);
Bm2 = [Bm2 Am2];
Am2 = [Am2; Am2];
Ac2 = Matrix(Am2);
Bc2 = Matrix(Bm2);

assert(all((Ac1*Bc1) == (Am1*Bm1)), 'Matrix Multiplication FAILED.');
assert(all((Ac2*Bc2) == (Am2*Bm2)), 'Matrix Multiplication FAILED.');
display('Matrix Multiplication successful.')

% ------------------------------------------------

Am = randi([-15 15], 6);
Ac = Matrix(Am);

Bm = randi([-15 15], [6,3]);
Bc = Matrix(Bm);

assert(all(Ac.' == Am.'), 'Matrix Transposition FAILED.');
assert(all(Bc.' == Bm.'), 'Matrix Transposition FAILED.');
display('Matrix Transposition successful.')

% % ------------------------------------------------

% Am = randi([1 15], 5);
% Ac = Matrix(Am);

% [Lm Um Pm] = lu(Am)
% [Lc Uc Pc] = Ac.LU_factor()

% assert(all(Lc == Lm), 'LU Decomposition FAILED.');
% assert(all(Uc == Um), 'LU Decomposition FAILED.');
% assert(all(Pc == Pm), 'LU Decomposition FAILED.');
% display('Lower-Upper Decomposition successful.')

% ------------------------------------------------

Am = randi([-15 15], 6);
Ac = Matrix(Am);

Bm = randi([-15 15], 3);
Bc = Matrix(Bm);


assert(all(inv(Ac) == inv(Am)), 'Matrix Inversion FAILED.');
assert(all(inv(Bc) == inv(Bm)), 'Matrix Inversion FAILED.');
display('Matrix Inversion successful.')

% ------------------------------------------------

Am = randi([1 15], 6);
Ac = Matrix(Am);

Bm = randi([1 15], 3);
Bc = Matrix(Bm);

dAm = det(Am); dAc = det(Ac);
dBm = det(Bm); dBc = det(Bc);

assert(abs(dAc - dAm) < 0.001, 'Determinant FAILED.');
assert(abs(dBc - dBm) < 0.001, 'Determinant FAILED.');
display('Calculating the Determinant successful.')

% ------------------------------------------------

Am1 = randi([-15 15], 6);
Ac1 = Matrix(Am1);
Am2 = randi([-15 15], 3);
Ac2 = Matrix(Am2);

Bm1 = randi([-15 15], 6);
Bc1 = Matrix(Bm1);
Bm2 = randi([-15 15], 3);
Bm2 = [Bm2 Bm2(:,2)];
Bc2 = Matrix(Bm2);

xm = linsolve(Am1, Bm1);
xc = ge_linsolve(Ac1, Bc1);

ym = linsolve(Am2, Bm2);
yc = ge_linsolve(Ac2, Bc2);

assert(all(xc == xm), 'Solving Linear System (GE) FAILED.');
assert(all(yc == ym), 'Solving Linear System (GE) FAILED.');
display('Solving Linear System (Gaussian Elimination) successful.')

% ------------------------------------------------

Am1 = randi([1 15], 3);
Ac1 = Matrix(Am1);

[evecs evals] = eig(Am1);
[dom_evec dom_eval] = eig_dominant(Ac1, 10);
assert(abs(dom_eval - max(max(evals))) < 0.0001, 'Dominant Eigenvalue FAILED');
display('Calculating Dominant Eigenvalue (Power Method) successful.')


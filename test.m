% Filename: test.m
% Description: Testing environment to validate custom operation implementations against the corresponding built-in MATLAB functions.


Am = [2 2 2; 3 3 3; 4 4 4];
Ac = Matrix(Am);

[Lm Um Pm] = lu(Am);
[Lc Uc Pc] = Ac.LU_factor();

Bm = magic(5);
Bc = Matrix(Bm);
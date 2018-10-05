
%************************************************
%                                               * 
%    Construct the Liouvillian superoperator    *
%        terms for the Heisenberg system.       *
%                                               *
%     (c) Dieter Jaksch and Stephen Clark       *
%                04.01.2012                     *
%                                               *
%************************************************

% Includes the Heisenberg Hamiltonian terms 
% within the coherent part of the Liouvillian.


function [system] = spin_liouvillian(n)

% Pauli matrices
ID = eye(2);
X = [0 1; 1 0];
Y = [0 -1i;1i 0];
Z = [1 0;0 -1];
Sp = [0 1; 0 0];
Sm = [0 0; 1 0];

coherent = 0; % Type flag for coherent Liouvillian terms.
incoherent = 1; % Type flag for incoherent Liouvillian terms.

% --- ONE SITE COHERENT TERMS ---

% - Z - for magnetic field
system(1).name = 'Z';
system(1).G = superoperators(Z,n,1,coherent);
system(1).type = 0; % One-site term.

% - Z - for tilted or random field (noted as electric field in setup file)
system(2).name = 'Zr';
system(2).G = superoperators(Z,n,1,coherent);
system(2).type = 0; % One-site term.

% --- TWO SITE COHERENT TERMS ---

% - XX -
system(3).name = 'XX'; 
system(3).G = superoperators(kron(X,X),n,2,coherent);
system(3).type = 1; % Two-site term.

% - YY -
system(4).name = 'YY';
system(4).G = superoperators(kron(Y,Y),n,2,coherent);
system(4).type = 1; % Two-site term.

% - ZZ -
system(5).name = 'ZZ';
system(5).G = superoperators(kron(Z,Z),n,2,coherent);
system(5).type = 1; % Two-site term.

% --- ONE SITE INCOHERENT TERMS ---

% - Z dephasing -
system(6).name = 'Z*';
system(6).G = superoperators(Z,n,1,incoherent);
system(6).type = 0; % One-site term.

% - S+ gain -
system(7).name = 'S+*';
system(7).G = superoperators(Sp,n,1,incoherent);
system(7).type = 0; % One-site term.

% - S- decay -
system(8).name = 'S-*';
system(8).G = superoperators(Sm,n,1,incoherent);
system(8).type = 0; % One-site term.
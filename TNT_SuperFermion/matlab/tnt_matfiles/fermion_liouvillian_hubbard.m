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


function [system] = fermion_liouvillian_hubbard(n)

% Pauli matrices
ID = eye(4);
nup = diag([0,1,0,1]);
ndn = diag([0,0,1,1]);
cup = [0,1,0,0;
       0,0,0,0;
       0,0,0,1;
       0,0,0,0];
cdup = cup';
cdn = [0,0,1,0;
       0,0,0,-1;
       0,0,0,0;
       0,0,0,0];
cddn = cdn';
sz = nup - ndn;
P = (-1)^(nup+ndn);

coherent = 0; % Type flag for coherent Liouvillian terms.
incoherent = 1; % Type flag for incoherent Liouvillian terms.

% --- ONE SITE COHERENT TERMS ---

% - U - Coulomb repulsion
system(1).name = 'Int';
system(1).G = superoperators(nup*ndn,n,2,coherent);
system(1).type = 0; % One-site term.

% --- TWO SITE COHERENT TERMS ---

% - UP_RL -
system(2).name = 'UP_RL'; 
system(2).G = superoperators(kron(cdup,ID)*kron(P,cup),n,4,coherent);
system(2).type = 1; % Two-site term.

% - UP_LR -
system(3).name = 'UP_LR'; 
system(3).G = superoperators(kron(P,cdup)*kron(cup,ID),n,4,coherent);
system(3).type = 1; % Two-site term.

% - DOWN_RL -
system(4).name = 'DOWN_RL'; 
system(4).G = superoperators(kron(cddn,ID)*kron(P,cdn),n,4,coherent);
system(4).type = 1; % Two-site term.

% - DOWN_LR -
system(5).name = 'DOWN_LR'; 
system(5).G = superoperators(kron(P,cddn)*kron(cdn,ID),n,4,coherent);
system(5).type = 1; % Two-site term.

% - ZZ -
system(6).name = 'ZZ'; 
system(6).G = superoperators(kron(sz,ID)*kron(ID,sz),n,4,coherent);
system(6).type = 1; % Two-site term.

% --- ONE SITE INCOHERENT TERMS ---

% - Gain_UP -
system(7).name = 'Gain_UP';
system(7).G = superoperators(cdup*(ID-ndn),n,2,incoherent);
system(7).type = 0; % One-site term.

% - LOSS_UP -
system(8).name = 'LOSS_UP';
system(8).G = superoperators(cup*(ndn),n,2,incoherent);
system(8).type = 0; % One-site term.

% - Gain_DOWN -
system(9).name = 'Gain_DOWN';
system(9).G = superoperators(cddn*(ID-nup),n,2,incoherent);
system(9).type = 0; % One-site term.

% - LOSS_UP -
system(10).name = 'LOSS_DOWN';
system(10).G = superoperators(cdn*(nup),n,2,incoherent);
system(10).type = 0; % One-site term.


%************************************************
%                                               * 
%    Constructs all the Hamiltonian matrices    *
%        required for Heisenberg system.        *
%                                               *
%     (c) Dieter Jaksch and Stephen Clark       *
%                23.02.2012                     *
%                                               *
%************************************************

% All terms required for the Heisenberg model are constructed here.

function [Sp,Sm,X,Y,Z,ID] = spin_hamiltonian(n)

% --- ONE SITE TERMS ---

S = (n - 1)/2; % The spin of each site.

ID = eye(n);
% Construct spin ladder operators :
Sm = zeros(n);
for m=-S:(S-1) 
  Sm(m+S+2,m+S+1) = sqrt(S*(S+1) - m*(m+1));
end;
Sp = zeros(n);
for m=-(S-1):S 
  Sp(m+S,m+S+1) = sqrt(S*(S+1) - m*(m-1));
end;
% Construct the spin X,Y and Z operators :
X = Sp + Sm;
Y = (Sp - Sm)/1i;
Z = diag([S:-1:-S])/S;


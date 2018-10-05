function Uf = fermionic_swap_gate
% Generates a fermionic swap gate for a Hubbard system:
% |1> = |vac>, |2> = |up>, |3> = |down> and |4> = |up,down>.

d = 4;
P = [1;-1;-1;1]; % Parities of the states on each site.

Uf = zeros(d^2,d^2);
for n=1:d
  for m=1:d
    in = (n-1)*d + m;
    out = (m-1)*d + n;
    % If the states on both sites are single-fermion states then add a
    % phase factor to account for anticommutivity.
    if P(n) == -1 && P(m) == -1
      % Add a (-1) for swapping states like |up>|dn> -> |dn>|up>.
      Uf(out,in) = -1; 
    else
      % Leave as is for swaps like |0>|up> -> |up>|0>.
      Uf(out,in) = 1;
    end;
  end;
end;
% -------------------------------------------
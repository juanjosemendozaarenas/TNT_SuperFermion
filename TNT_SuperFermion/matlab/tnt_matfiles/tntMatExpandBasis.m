%======================================================================
%> @ingroup matscripts
%> Expands a matrix representing an operator for a single species into a matrix representing a multi-species operator
%>
%> Returns the result in a cell having a number of elements equal to the number of species.
%> If the original operator was \f$d\times\f$d, the new matrices will be \f$d^N\timesd^N\f$.
%>
%> @param ssop Single species operator
%> @param N The number of species
%>
%> @retval exop Cell array containing the operator for each species
%======================================================================
function ex_op = tntMatExpandBasis(ssop, N)

d = size(ssop,1);

ex_op = cell(N,1);

for n = 1:N
    ex_op{n} = kron(kron(eye(d^(n-1)),ssop),eye(d^(N-n)));
end

end
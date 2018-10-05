%======================================================================
%> @ingroup matscripts
%> Creates the tntSystem structure and initialises it with the basis operator structure, where the basis operator has been provided as a matrix.
%> If there is a global physical symmetry (e.g. total number of bosons are being conserved), then the associated quantum numbers should be provided. Otherwise simply pass an empty array.
%> In this case the basis operator should be diagonal
%>
%> @param op Matrix representing the basis operator
%> @param qnums Quantum numbers for a physical leg of the operator
%>
%> @retval sx A struct
%======================================================================

function tntSystem = tntMatCreateBasisOp(op,qnums)
    d = size(op,1);
if (~isempty(qnums))
    tensor.elems_type = 'blocks';
    if (1 == size(qnums,1))
        [qnumssorted,indmap] = sort(qnums);  %qnumssorted(j) = qnums(i)  j = indmap(i)
        % find the number of blocks
        numblocks = 1;
        indblock = 0;
        for j=1:d
            if ((j > 1) && qnumssorted(j) ~= qnumssorted(j-1))
                numblocks = numblocks + 1;
                indblock = 1;
            else
                indblock = indblock+1;
            end
            tensor.elems.qn_tot(numblocks) = qnumssorted(j);
            i = find(indmap==j);
            tensor.elems.vals{numblocks}(indblock,indblock) = op(i,i);
            tensor.elems.indmapr(1,i) = numblocks-1;
            tensor.elems.indmapr(2,i) = indblock-1;
        end
        tensor.elems.indmapc = tensor.elems.indmapr;
    else
        tensor.elems.qn_tot = qnums;
        tensor.elems.vals = num2cell(diag(op));
        tensor.elems.indmapr = [0:(d -1);zeros(1,d)];
        tensor.elems.indmapc = [0:(d-1);zeros(1,d)];
    end
    tensor.elems.rowlegs = 0;
    tensor.elems.collegs = 1;
    tensor.qn_info.qn_dir = [1,-1];
    tensor.qn_info.qn_index = {qnums,qnums};
else
    tensor.elems_type = 'values';
    tensor.elems = op;
    tensor.qn_info.qn_dir = [0,0];
    tensor.qn_info.qn_index = {[],[]};
end

tensor.dims = size(op);

tntSystem.basisOp.tensor = tensor;
tntSystem.basisOp.ids = 'DU';
tntSystem.basisOp.indices = {0,1};

end



    



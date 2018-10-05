

%======================================================================
%> @ingroup matscripts
%> Creates a structure that represents a product MPS network, that can then be loaded uses tntNetworksLoad().
%> If the passed cfg is not normalized, here it is left like that
%>
%> @param cfg A cell, each entry containing a vector representing the state on each site
%> @param qnums The quantum numbers for the physical leg, or an empty array if there is no global symmetry
%>
%> @retval wf A structure representing the product MPS
%======================================================================

function wf = tntMatCreateProdMpsNotNormalized(cfg,qnums)

d = length(cfg{1});
L = length(cfg);
qcurr = zeros(size(qnums,1),1);

for loop = 1:L
    if (~isempty(qnums))
        ind = find(cfg{loop},1);
        tensor.elems_type = 'blocks';
        tensor.elems.qn_tot = qcurr+qnums(:,ind);
        tensor.elems.vals{1} = 1;
        tensor.elems.rowlegs = [0,1];
        tensor.elems.collegs = 2;
        tensor.elems.indmapr = [-ones(1,d);zeros(1,d)];
        tensor.elems.indmapr(1,ind) = 0;
        tensor.elems.indmapc = [0;0];
        tensor.qn_info.qn_dir = [1,1,-1];
        tensor.qn_info.qn_index = {qcurr,qnums,qcurr+qnums(:,ind)};
        qcurr = qcurr+qnums(:,ind);
    else
        tensor.elems_type = 'values';
        tensor.elems = cfg{loop};
        tensor.qn_info.qn_dir = [0,0,0];
        tensor.qn_info.qn_index = {[],[],[]};
    end
    tensor.dims = [1,d,1];
    
    wf.nodes(loop).tensor = tensor;
    wf.nodes(loop).ids = 'LRD';
    wf.nodes(loop).indices = {0,2,1};
end

wf.start = 0;
wf.start_leg = 'L';
wf.end = L-1;
wf.end_leg = 'R';

wf.connections = cell(L,L);

for loop=2:L
    wf.connections{loop-1,loop} = 'RL';
end

end
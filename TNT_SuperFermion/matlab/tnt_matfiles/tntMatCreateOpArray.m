%======================================================================
%> @ingroup matscripts
%> Turns a matrix into a node structure representing an operator in the MPS library.
%> To create a single node, call the function passing a cell of length 1 containing the matrix. The node can then be loaded using tntNodesLoad().
%>
%> To create a node array, call the function passing a cell array, each entry contraining a matrix. The node can then be loaded using tntNodeArraysLoad().
%>
%> The nodes are created have two legs according to the labelling of the MPS library i.e. columns correspond to the upwards facing leg and are labelled "U", rows correspond to the downwards facing leg and are labelled "D".
%>
%> @param ops A cell array contraining matrices for the operators.
%>
%> @retval narr A structure representing a singe node or a node array.
%======================================================================
function narr = tntMatCreateOpArray(ops)

if (isempty(ops))
    narr = struct('tensor',{},'ids',{},'indices',{});
end

for loop = 1:length(ops)
    tensor.elems_type = 'values';
    tensor.elems = ops{loop};
    tensor.dims = size(ops{loop});
    tensor.qn_info.qn_dir = [0,0];
    tensor.qn_info.qn_index = {[],[]};
    narr(loop).tensor = tensor;
    narr(loop).ids = 'DU';
    narr(loop).indices = {0,1};
end

end
%======================================================================
%> @ingroup matscripts
%> Turns a cell array of operators into a functional node representing a
%> sinlge site MPS operator. Once loaded into the TNT library, parameters on
%> the node can be set to change the final value.
%> 
%> @param ops A cell array of the matrix operators.
%>
%> @param func A string representing the function for the node - either
%> 'sum' or 'exp' are currently supported.
%>
%> @retval node A structure representing the functional node.
%======================================================================
function node = tntMatCreateFunctionalOp(ops,func)

tensor.elems_type = 'params';
tensor.elems.function = func;
tensor.elems.operators = ops;
tensor.dims = size(ops{1});
tensor.qn_info.qn_dir = [0,0];
tensor.qn_info.qn_index = {[],[]};
node.tensor = tensor;
node.ids = 'DU';
node.indices = {0,1};

end
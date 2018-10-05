% This function is for creating two-site gates with the structure in the
% new TNT library. In this case we have four legs (two up and two down).
function narr = tntMatCreate2SitesOpArray(ops)

for loop = 1:length(ops)
    tensor.elems_type = 'values';
    tensor.elems = ops{loop};
    tensor.dims = size(ops{loop});
    tensor.qn_info.qn_dir = [0,0,0,0];
    tensor.qn_info.qn_index = {[],[],[],[]};
    narr(loop).tensor = tensor;
    narr(loop).ids = 'DEUV';
    narr(loop).indices = {0,1,2,3};
end

end
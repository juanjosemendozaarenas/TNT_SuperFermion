% This function is for creating tensors with a single leg 

function narr = tntMatCreateOneLegOpArray(ops)

for loop = 1:length(ops)
    tensor.elems_type = 'values';
    tensor.elems = ops{loop};
    tensor.dims = size(ops{loop},2); % To remove the information about the leg with dimension 1
    tensor.qn_info.qn_dir = [0];
    tensor.qn_info.qn_index = {[]};
    narr(loop).tensor = tensor;
    narr(loop).ids = 'U';
    narr(loop).indices = {0};
end

end
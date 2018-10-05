% This function is for creating tensors with three legs 

function narr = tntMatCreateThreeLegsOpArray(ops)

for loop = 1:length(ops)
    tensor.elems_type = 'values';
    tensor.elems = ops{loop};
    tensor.dims = size(ops{loop});
    tensor.qn_info.qn_dir = [0,0,0];
    tensor.qn_info.qn_index = {[],[],[]};
    narr(loop).tensor = tensor;
    narr(loop).ids = 'LRU';
    narr(loop).indices = {0,1,2}; % This has to agree with the order of the dimensions!!
end

end
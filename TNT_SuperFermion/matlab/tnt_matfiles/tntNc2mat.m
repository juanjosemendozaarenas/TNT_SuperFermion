function tntNc2mat(ncname, matoutput)
%#ok<*NASGU> 
%#ok<*ASGLU>
ncid = netcdf.open(ncname,'NOWRITE');

% ----- Reading the system information ----- %
verid = netcdf.inqVarID(ncid,'tntLibType');
tntLibType = netcdf.getVar(ncid,verid);
tntLibVersion = tntgetAtt(ncid,verid,'tntLibVersion');
gid = netcdf.inqNcid(ncid,'system');
verid = netcdf.inqVarID(gid,'sysinfo');
tntSystem_nc = struct;
tntSystem_nc.sysnum = tntgetAtt(gid,verid,'sysnum');
tntSystem_nc.symm_type = tntgetAtt(gid,verid,'symm_type');
tntSystem_nc.symm_num_qn = tntgetAtt(gid,verid,'symm_num_qn');
tntSystem_nc.zero_tol = tntgetAtt(gid,verid,'zero_tol');
tntSystem_nc.abs_trunc_tol = tntgetAtt(gid,verid,'abs_trunc_tol');
tntSystem_nc.rel_trunc_tol = tntgetAtt(gid,verid,'rel_trunc_tol');
tntSystem_nc.svdtype = tntgetAtt(gid,verid,'svdtype');
tntSystem_nc.maxeigiter = tntgetAtt(gid,verid,'maxeigiter');
tntSystem_nc.reshape_reuse = tntgetAtt(gid,verid,'reshape_reuse');
tntSystem_nc.trunc_err_func = tntgetAtt(gid,verid,'trunc_err_func');
basisopid = netcdf.inqGrps(gid);
if (~isempty(basisopid))
    tntSystem_nc.basisOp = tntNcGetNode(basisopid);
else
    tntSystem_nc.basisOp = [];
end
tntSystem = tntSystem_nc; 
tntFileInfo = ['NetCDF file ',ncname,' converted to MATLAB format']; 
save(matoutput, 'tntSystem', 'tntLibType', 'tntLibVersion','tntFileInfo');

% ----- Reading parameters in the output file ------
gid = netcdf.inqNcid(ncid,'parameters');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    [varname, vardata] = tntNcGetData(varids(j)); 
    eval([varname,' = transpose(vardata);']);
    save(matoutput,varname,'-append');
    disp(['Converted parameter ',varname]);
end

% ----- Reading arrays in the output file ------
gid = netcdf.inqNcid(ncid,'arrays');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    [varname, vardata] = tntNcGetData(varids(j));
    eval([varname,' = vardata;']);
    save(matoutput, varname, '-append');
    disp(['Converted array ',varname]);
end

% ----- Reading strings in the output file ------
gid = netcdf.inqNcid(ncid,'strings');
varids = netcdf.inqVarIDs(gid);
for j = 1:length(varids)
    varname = netcdf.inqVar(gid,varids(j));
    vardata = tntgetAtt(gid,varids(j),'value');
    eval([varname,' = vardata;']);
    save(matoutput,varname,'-append');
    disp(['Converted string ',varname]);
end

% ----- Reading nodes in the output file ------
gid = netcdf.inqNcid(ncid,'nodes');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    varname = netcdf.inqGrpName(varids(j));
    nodestruct = tntNcGetNode(varids(j));
    eval([varname,' = nodestruct;']);
    save(matoutput, varname, '-append');
    disp(['Converted node ',varname]);
end

% ----- Reading nodes arrays in the output file ------
gid = netcdf.inqNcid(ncid,'node arrays');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    varname = netcdf.inqGrpName(varids(j));
    nodeids = netcdf.inqGrps(varids(j));
    nodestruct = tntNcGetNode(nodeids);
    eval([varname,' = nodestruct;']);
    save(matoutput, varname, '-append');
    disp(['Converted node array ',varname]);
end

%----- Reading expectation value operators in the output file ------
gid = netcdf.inqNcid(ncid,'exops');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    varname = netcdf.inqGrpName(varids(j));
    exop = struct;
    
    gid = netcdf.inqNcid(varids(j),'os_operators');
    ops = netcdf.inqGrps(gid);
    exop.os_operators = tntNcGetNode(ops);
    exop.os_labels = cell(1, length(ops));
    for o = 1:length(ops)
       exop.os_labels{o} = netcdf.inqGrpName(ops(o));
    end
       
    gid = netcdf.inqNcid(varids(j),'nn_operators');
    ops = netcdf.inqGrps(gid);
    exop.nn_operators = tntNcGetNode(ops);
    exop.nn_labels = cell(1, length(ops)/2);
    for o = 2:2:length(ops)
       exop.nn_labels{o/2} = netcdf.inqGrpName(ops(o));
       exop.nn_labels{o/2} = exop.nn_labels{o/2}(1:end-1);
    end
    
    gid = netcdf.inqNcid(varids(j),'cs_operators');
    ops = netcdf.inqGrps(gid);
    exop.cs_operators = tntNcGetNode(ops);
    exop.cs_labels = cell(1, length(ops)/2);
    for o = 2:2:length(ops)
       exop.cs_labels{o/2} = netcdf.inqGrpName(ops(o));
       exop.cs_labels{o/2} = exop.cs_labels{o/2}(1:end-1);
    end
    
    gid = netcdf.inqNcid(varids(j),'ap_operators');
    ops = netcdf.inqGrps(gid);
    exop.ap_operators = tntNcGetNode(ops);
    exop.ap_labels = cell(1, length(ops)/2);
    for o = 2:2:length(ops)
       exop.ap_labels{o/2} = netcdf.inqGrpName(ops(o));
       exop.ap_labels{o/2} = exop.ap_labels{o/2}(1:end-1);
    end

    eval([varname,' = exop;']);
    save(matoutput, varname, '-append');
    disp(['Converted exop ',varname]);
end

% ----- Reading networks in the output file ------
gid = netcdf.inqNcid(ncid,'networks');
varids = netcdf.inqGrps(gid);
for j = 1:length(varids)
    varname = netcdf.inqGrpName(varids(j));
    nodeids = netcdf.inqGrps(netcdf.inqNcid(varids(j),'nodes'));
    conngid = netcdf.inqNcid(varids(j),'connections');
    connid = netcdf.inqVarID(conngid,'conn_nodes');
    networkstruct.nodes = tntNcGetNode(nodeids);
    numnodes = length(networkstruct.nodes);
    dimid = netcdf.inqDimID(conngid,'connecting legs');
    [~, dimlen] = netcdf.inqDim(conngid,dimid);
    networkstruct.connections = cell(numnodes);
    for m = 1:numnodes
        for n = m:numnodes
            str = netcdf.getVar(conngid, connid,[0,m-1,n-1],[dimlen,1,1]);
            networkstruct.connections{m,n} = str(isstrprop(str, 'alphanum'))'; 
        end
    end
    networkstruct.start = tntgetAtt(conngid, connid, 'start');
    networkstruct.start_leg = tntgetAtt(conngid, connid, 'start_leg');
    networkstruct.end = tntgetAtt(conngid, connid, 'end');
    networkstruct.end_leg = tntgetAtt(conngid, connid, 'end_leg');
    eval([varname,' = networkstruct;']);
    save(matoutput, varname, '-append');
    disp(['Converted network ',varname]);
end

netcdf.close(ncid);
end

function [varname,vardata] = tntNcGetData(varid)
    varname = netcdf.inqGrpName(varid);
    vardata = netcdf.getVar(varid,netcdf.inqVarID(varid,'real_part'));
    if (length(netcdf.inqVarIDs(varid)) > 1) 
        vardata = vardata + 1i*netcdf.getVar(varid,netcdf.inqVarID(varid,'imaginary_part')); 
    %elseif any(isinteger(vardata))
    %    vardata = double(vardata);
    end
end

function nodestruct = tntNcGetNode(nodeids)

nodestruct = struct;

for k = 1:length(nodeids)
        leginfo = netcdf.inqVarID(nodeids(k),'legs');
        nodestruct(k).ids = tntgetAtt(nodeids(k),leginfo,'ids');
        
        numlegs = netcdf.getVar(nodeids(k),netcdf.inqVarID(nodeids(k),'legs'));
        nodestruct(k).indices = cell(1,numlegs);
        for l=1:numlegs
            nodestruct(k).indices{l} = tntgetAtt(nodeids(k),leginfo,['leg_',num2str(l-1),'_indices']);
        end
        
        tensorid = netcdf.inqNcid(nodeids(k),'tensor');
        tensorinfo = netcdf.inqVarID(tensorid,'info');
        nodestruct(k).tensor.elems_type = tntgetAtt(tensorid,tensorinfo,'elems_type');
        nodestruct(k).tensor.dims = tntgetAtt(tensorid,tensorinfo,'dims');
        nodestruct(k).tensor.qn_info.qn_dir = tntgetAtt(tensorid,tensorinfo,'qn_dir');
        nodestruct(k).tensor.qn_info.qn_index = cell(1,numlegs);
        qnsymmnum = tntgetAtt(tensorid,tensorinfo,'qn_symmnum');
        if (qnsymmnum && any(nodestruct(k).tensor.qn_info.qn_dir))
            for i=1:length(nodestruct(k).tensor.qn_info.qn_dir)
                    nodestruct(k).tensor.qn_info.qn_index{i} = tntgetAtt(tensorid,tensorinfo,['qn_index_',num2str(i-1)]);
                    nodestruct(k).tensor.qn_info.qn_index{i} = reshape(nodestruct(k).tensor.qn_info.qn_index{i}, [qnsymmnum, nodestruct(k).tensor.dims(i)]);
            end
        end
        elemsid = netcdf.inqNcid(tensorid,'elems');
        
        if (strcmp(nodestruct(k).tensor.elems_type,'values'))
            elems.vals = tntNcGetTensorValues(elemsid, 'vals');
        elseif (strcmp(nodestruct(k).tensor.elems_type,'blocks'))
            infoid = netcdf.inqVarID(elemsid,'info');
            numblocks = netcdf.getVar(elemsid,infoid);
            elems.vals = cell(1,numblocks);
            for b = 1:numblocks
                elems.vals{b} = tntNcGetTensorValues(elemsid, ['block_' num2str(b-1)]);
            end
            elems.qn_tot = reshape(tntgetAtt(elemsid,infoid,'qn_tot'),[qnsymmnum,numblocks]);
            elems.rowlegs = tntgetAtt(elemsid,infoid,'rowlegs');
            elems.collegs = tntgetAtt(elemsid,infoid,'collegs');
            elems.indmapr = tntgetAtt(elemsid,infoid,'indmapr');
            elems.indmapc = tntgetAtt(elemsid,infoid,'indmapc');
            elems.indmapr = reshape(elems.indmapr,[2,length(elems.indmapr)/2]);
            elems.indmapc = reshape(elems.indmapc,[2,length(elems.indmapc)/2]);
        end
        nodestruct(k).tensor.elems = elems;
end

end

function elemvals = tntNcGetTensorValues(elemsid, varname)
    
valsid = netcdf.inqVarID(elemsid,varname);
rowsize = tntgetAtt(elemsid,valsid,'rowsize');
colsize = tntgetAtt(elemsid,valsid,'colsize');
comp = tntgetAtt(elemsid,valsid,'comp');
elemvals = netcdf.getVar(elemsid,valsid,[0 0 0],[comp+1 rowsize colsize]);
if (comp)
    elemvals = elemvals(1,:,:) + 1i*elemvals(2,:,:);
end
elemvals = reshape(elemvals,[rowsize colsize]);

end

function convertedval = tntgetAtt(gid,varid,varname)
    convertedval = netcdf.getAtt(gid,varid,varname);
    
    if (any(isinteger(convertedval)))
        convertedval = double(convertedval);
    end
end
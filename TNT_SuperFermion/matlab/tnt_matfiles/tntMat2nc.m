function tntMat2nc(matname, ncoutput)
%#ok<*NASGU> 
%#ok<*ASGL

if (exist(ncoutput,'file'))
    delete(ncoutput)
end

disp(['-- Creating file ' ncoutput ' --']);

% Create the file and all the groups needed;
ncid = netcdf.create(ncoutput,'NETCDF4');

systemid = netcdf.defGrp(ncid,'system');
arraysid = netcdf.defGrp(ncid,'arrays');
parametersid = netcdf.defGrp(ncid,'parameters');
nodesid = netcdf.defGrp(ncid,'nodes');
nodearrsid = netcdf.defGrp(ncid,'node arrays');
networksid = netcdf.defGrp(ncid,'networks');
exopsid = netcdf.defGrp(ncid,'exops');
stringsid = netcdf.defGrp(ncid,'strings');
% define the dimension for all variables
stringdimid = netcdf.defDim(stringsid,'string length',1);

% Find all the variables in the matlab file
mvars = who('-file',matname);
load(matname);

% If there is no library type already set (i.e. if this is from a manually
% created init file), set the librry type to -1, then save this variable to
% the file
if (~exist('tntLibType','var'))
    tntLibType = -1;
end
if (~exist('tntLibVersion','var'))
    tntLibVersion = 'manually generated file';
end
verid = netcdf.defVar(ncid,'tntLibType','NC_INT',[]);
netcdf.putVar(ncid,verid,tntLibType);
netcdf.putAtt(ncid,verid,'tntLibVersion',tntLibVersion);

% Go through each variable
for loop=1:length(mvars)
    mvar = eval(mvars{loop});
    % determine what the variable is
    if (strcmp('tntLibType',mvars{loop}) || strcmp('tntLibVersion',mvars{loop}))
        %do nothing
    elseif (strcmp('tntSystem',mvars{loop}))
%% SYSTEM INFORMATION
        disp('Saving system information to netcdf file');
        verid = netcdf.defVar(systemid,'sysinfo','NC_INT',[]);
        dummyvar = 0;
        netcdf.putVar(systemid,verid,dummyvar);
        netcdf.putAtt(systemid,verid,'sysnum',tntSystem.sysnum);
        netcdf.putAtt(systemid,verid,'symm_type',tntSystem.symm_type);
        netcdf.putAtt(systemid,verid,'symm_num_qn',tntSystem.symm_num_qn);
        netcdf.putAtt(systemid,verid,'zero_tol',tntSystem.zero_tol);
        netcdf.putAtt(systemid,verid,'abs_trunc_tol',tntSystem.abs_trunc_tol);
        netcdf.putAtt(systemid,verid,'rel_trunc_tol',tntSystem.rel_trunc_tol);
        netcdf.putAtt(systemid,verid,'svdtype',tntSystem.svdtype);
        netcdf.putAtt(systemid,verid,'maxeigiter',tntSystem.maxeigiter);
        netcdf.putAtt(systemid,verid,'reshape_reuse',tntSystem.reshape_reuse);
        netcdf.putAtt(systemid,verid,'trunc_err_func', tntSystem.trunc_err_func);
        if (~isempty(tntSystem.basisOp))
            basisopid = netcdf.defGrp(systemid,'basisOp');
            tntNcWriteNode(basisopid,tntSystem.basisOp);
        end

    elseif (ischar(mvar))
        disp(['Saving string variable ',mvars{loop},' to netcdf file']);
%% STRING VARIABLES
        % Create a variable for the string 
        varid = netcdf.defVar(stringsid,mvars{loop},'NC_UINT',stringdimid);
        netcdf.putVar(stringsid,varid,length(mvar));
        netcdf.putAtt(stringsid,varid,'value',mvar);

    elseif (isnumeric(mvar))
        if (length(mvar) == 1)
%% PARAMETERS
            % Create a group for the parameter variable 
            gid = netcdf.defGrp(parametersid,mvars{loop});
            dimid = [];
            dimid(1) = netcdf.defDim(gid,'value',1);
            dimid(2) = netcdf.defDim(gid,'update',netcdf.getConstant('NC_UNLIMITED'));
         
            if (isreal(mvar))
                disp(['Saving real parameter ',mvars{loop},' to netcdf file']);
                if (isinteger(mvar))
                    varid = netcdf.defVar(gid,'real_part','NC_INT',dimid);
                else
                    varid = netcdf.defVar(gid,'real_part','NC_DOUBLE',dimid);
                end
                netcdf.putVar(gid,varid,[0 0],[1 1],mvar);
            else 
                disp(['Saving complex parameter ',mvars{loop},' to netcdf file']);
                varid = netcdf.defVar(gid,'real_part','NC_DOUBLE',dimid);
                netcdf.putVar(gid,varid,[0 0],[1 1],real(mvar));
                varid = netcdf.defVar(gid,'imaginary_part','NC_DOUBLE',dimid);
                netcdf.putVar(gid,varid,[0 0],[1 1],imag(mvar));
            end   
        else
%% ARRAYS
            % Create a group for the array variable 
            gid = netcdf.defGrp(arraysid,mvars{loop});
            rowsize = size(mvar,1);
            colsize = size(mvar,2);
            usize = numel(mvar)/(rowsize*colsize);
            dimid = [];
            dimid(1) = netcdf.defDim(gid,'rows',rowsize);
            dimid(2) = netcdf.defDim(gid,'columns',colsize);
            dimid(3) = netcdf.defDim(gid,'update',netcdf.getConstant('NC_UNLIMITED'));
            
            if (~isempty(mvar))
                if (isreal(mvar))
                    disp(['Saving real array ',mvars{loop},' to netcdf file']);
                    if (isinteger(mvar))
                        varid = netcdf.defVar(gid,'real_part','NC_INT',dimid);
                    else
                        varid = netcdf.defVar(gid,'real_part','NC_DOUBLE',dimid);
                    end
                    if numel(mvar); netcdf.putVar(gid,varid,[0 0 0],[rowsize colsize usize],mvar); end;
                else 
                    varid = netcdf.defVar(gid,'real_part','NC_DOUBLE',dimid);
                    if numel(mvar); netcdf.putVar(gid,varid,[0 0 0], [rowsize colsize usize], real(mvar)); end;
                    varid = netcdf.defVar(gid,'imaginary_part','NC_DOUBLE',dimid);
                    if numel(mvar); netcdf.putVar(gid,varid,[0 0 0], [rowsize colsize usize], imag(mvar)); end;
                end  
            else 
                disp(['Saving empty array ',mvars{loop},' to netcdf file']);
            end
        end  
    elseif (isstruct(mvar))
        if (isfield(mvar,'tensor'))
            if (length(mvar) == 1)
%% NODES
                disp(['Saving node ',mvars{loop},' to netcdf file']);
                % Create a group for the node variable 
                gid = netcdf.defGrp(nodesid,mvars{loop});
                tntNcWriteNode(gid,mvar);
            else
%% NODE ARRAYS
                disp(['Saving node array ',mvars{loop},' to netcdf file']);
                % Create a group for the node array variable 
                arrid = netcdf.defGrp(nodearrsid,mvars{loop});
                for k=1:length(mvar)
                    gid = netcdf.defGrp(arrid,['node_',num2str(k-1)]);
                    tntNcWriteNode(gid,mvar(k));
                end
            end
        elseif (isfield(mvar,'connections'))
            
%% NETWORKS
            disp(['Saving network ',mvars{loop},' to netcdf file']);
            
            % Define a group for the variable 
            gid = netcdf.defGrp(networksid, mvars{loop});
            
            numnodes = length(mvar.nodes);
            numc = zeros(numnodes,numnodes);
            
            % Convert the nodes in the network
            nnodeids = netcdf.defGrp(gid,'nodes');
            for k = 1:numnodes
                nnodeid = netcdf.defGrp(nnodeids,['node_',num2str(k-1)]);
                tntNcWriteNode(nnodeid,mvar.nodes(k));
            end
            
            % Convert the connections in the network
            conngid = netcdf.defGrp(gid,'connections');
            dimid = [];
            dimid(3) = netcdf.defDim(conngid,'node1',numnodes);
            dimid(2) = netcdf.defDim(conngid,'node2',numnodes);
            dimid(1) = netcdf.defDim(conngid,'connecting legs',netcdf.getConstant('NC_UNLIMITED'));
            
            connid = netcdf.defVar(conngid,'conn_nodes','NC_CHAR',dimid);
            for m = 1:numnodes
                for n = m:numnodes
                    numc(n,m) = length(mvar.connections{m,n})/2;
                    numc(m,n) = numc(n,m);
                    if (numc(n,m))
                        strvar = mvar.connections{m,n};
                        netcdf.putVar(conngid, connid, [0 m-1 n-1], [2*numc(n,m) 1  1], strvar);
                    end
                end
            end
            netcdf.putAtt(conngid, connid, 'numc', int32(numc));
            netcdf.putAtt(conngid, connid, 'start', int32(mvar.start));
            netcdf.putAtt(conngid, connid, 'start_leg', mvar.start_leg);
            netcdf.putAtt(conngid, connid, 'end', int32(mvar.end));
            netcdf.putAtt(conngid, connid, 'end_leg', mvar.end_leg);
        elseif (isfield(mvar,'os_operators'))
%% EXPECTATION OPERATORS
            disp(['Saving expectation operators ',mvars{loop},' to netcdf file']);
            
            % Define a group for the variable
            varid = netcdf.defGrp(exopsid,mvars{loop});
    
            gid = netcdf.defGrp(varid,'os_operators');
            for o = 1:length(mvar.os_labels)
                opid = netcdf.defGrp(gid,mvar.os_labels{o});
                tntNcWriteNode(opid,mvar.os_operators(o));
            end

            gid = netcdf.defGrp(varid,'nn_operators');
            for o = 1:length(mvar.nn_labels)
                opid = netcdf.defGrp(gid,[mvar.nn_labels{o} 'L']);
                tntNcWriteNode(opid,mvar.nn_operators(2*o-1));
                
                opid = netcdf.defGrp(gid,[mvar.nn_labels{o} 'R']);
                tntNcWriteNode(opid,mvar.nn_operators(2*o));
            end

            gid = netcdf.defGrp(varid,'cs_operators');
            for o = 1:length(mvar.cs_labels)
                opid = netcdf.defGrp(gid,[mvar.cs_labels{o} 'L']);
                tntNcWriteNode(opid,mvar.cs_operators(2*o-1));
                
                opid = netcdf.defGrp(gid,[mvar.cs_labels{o} 'R']);
                tntNcWriteNode(opid,mvar.cs_operators(2*o));
            end
            
            gid = netcdf.defGrp(varid,'ap_operators');
            for o = 1:length(mvar.ap_labels)
                opid = netcdf.defGrp(gid,[mvar.ap_labels{o} 'L']);
                tntNcWriteNode(opid,mvar.ap_operators(2*o-1));
                
                opid = netcdf.defGrp(gid,[mvar.ap_labels{o} 'R']);
                tntNcWriteNode(opid,mvar.ap_operators(2*o));
            end
        else
            disp(['Unknown variable type for ',mvars{loop}]);
        end  
    else
        disp(['Unknown variable type for ',mvars{loop}]);
    end
    
end
netcdf.close(ncid);
end
    
function tntNcWriteNode(gid,nodestruct)

        % write info about node legs to the variable and attributes
        leginfo = netcdf.defVar(gid,'legs','NC_UINT',[]);
        numlegs = length(nodestruct.ids);
        netcdf.putVar(gid,leginfo,numlegs);
        netcdf.putAtt(gid,leginfo,'ids',nodestruct.ids); 
        for l=1:numlegs
            netcdf.putAtt(gid,leginfo,['leg_',num2str(l-1),'_indices'], uint32(nodestruct.indices{l}));
        end
        
        % write tensor info to a new group
        tensorid = netcdf.defGrp(gid,'tensor');
        tensorinfo = netcdf.defVar(tensorid,'info','NC_UINT',[]);
        netcdf.putAtt(tensorid,tensorinfo,'elems_type', nodestruct.tensor.elems_type);
        netcdf.putAtt(tensorid,tensorinfo,'dims', int32(nodestruct.tensor.dims));
        
        % write quantum number information as attributes
        netcdf.putAtt(tensorid,tensorinfo,'qn_dir', int32(nodestruct.tensor.qn_info.qn_dir));
        if all(nodestruct.tensor.qn_info.qn_dir == 0)
            qnsymmnum = 0;
        else
            qnsymmnum = size(nodestruct.tensor.qn_info.qn_index{find(nodestruct.tensor.qn_info.qn_dir,1)},1);
        end
        netcdf.putAtt(tensorid,tensorinfo,'qn_symmnum', int32(qnsymmnum));
        if (qnsymmnum)
            for i=1:length(nodestruct.tensor.qn_info.qn_dir)
                    if (nodestruct.tensor.qn_info.qn_dir(i))
                        %reshape(nodestruct(k).tensor.qn_info.qn_index{i}, [qnsymmnum, nodestruct(k).tensor.dims(i)]);
                        netcdf.putAtt(tensorid,tensorinfo,['qn_index_',num2str(i-1)], int32(nodestruct.tensor.qn_info.qn_index{i}));
                    end
            end
        end
        
        % Create a group for tensor elements
        elemsid = netcdf.defGrp(tensorid,'elems');
        dimid = [];
        dimid(1) = netcdf.defDim(elemsid,'complex',netcdf.getConstant('NC_UNLIMITED'));
        dimid(2) = netcdf.defDim(elemsid,'rows',netcdf.getConstant('NC_UNLIMITED'));
        dimid(3) = netcdf.defDim(elemsid,'columns',netcdf.getConstant('NC_UNLIMITED'));
        
        elems = nodestruct.tensor.elems;
        infoid = netcdf.defVar(elemsid,'info','NC_UINT',[]);
        
        if (strcmp(nodestruct.tensor.elems_type,'values'))
            
            netcdf.putVar(elemsid,infoid,1);
            if (isstruct(elems)) 
                tntNcWriteTensorValues(elemsid, dimid, 'vals', elems.vals);
            else %old format of init files
                tntNcWriteTensorValues(elemsid, dimid, 'vals', elems);
            end
        elseif (strcmp(nodestruct.tensor.elems_type,'blocks'))
            
            numblocks = length(elems.vals);
            netcdf.putVar(elemsid,infoid,numblocks);
            for b = 1:numblocks
                tntNcWriteTensorValues(elemsid, dimid, ['block_' num2str(b-1)], elems.vals{b});
            end
            % put information as attributes
            netcdf.putAtt(elemsid,infoid,'qn_tot',int32(elems.qn_tot(:)));
            netcdf.putAtt(elemsid,infoid,'rowlegs',uint32(elems.rowlegs));
            netcdf.putAtt(elemsid,infoid,'collegs',uint32(elems.collegs));
            netcdf.putAtt(elemsid,infoid,'indmapr',int32(elems.indmapr(:)));
            netcdf.putAtt(elemsid,infoid,'indmapc',int32(elems.indmapc(:)));

        end


end

function tntNcWriteTensorValues(elemsid, dimid, varname, elemvals)

valsid = netcdf.defVar(elemsid,varname,'NC_DOUBLE',dimid);

comp = ~isreal(elemvals);
rowsize = size(elemvals,1);
colsize = size(elemvals,2);

if (comp)
    var = zeros(rowsize,colsize,comp+1);
    var(:,:,1) = real(elemvals);
    var(:,:,2) = imag(elemvals);
    var = permute(var,[3,1,2]);
else
    var = reshape(elemvals,[1 rowsize, colsize]);
end

netcdf.putVar(elemsid,valsid,[0 0 0],[comp+1 rowsize colsize], var);
netcdf.putAtt(elemsid,valsid,'comp',int32(comp));
netcdf.putAtt(elemsid,valsid,'rowsize',int32(rowsize));
netcdf.putAtt(elemsid,valsid,'colsize',int32(colsize));
end

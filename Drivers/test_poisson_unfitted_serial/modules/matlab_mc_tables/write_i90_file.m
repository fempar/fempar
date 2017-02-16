function fc = write_i90_file(file_name, mc_ncases, mc_max_sub_cells, mc_num_sub_cells_per_case,...
    mc_subcells_per_case, mc_inout_subcells_per_case, mc_num_nodes_per_subcell,mc_num_cut_edges_per_case,elem_type)

% Falta el num d'intersection points

% Open file
fid = fopen(file_name,'w');


fprintf(fid,'! Look up tables for %s\n',elem_type);
fprintf(fid,'! This file has been automatically generated in Matlab.\n');
fprintf(fid,'! Do not modify this file by hand!\n');

fprintf(fid,'\n');

fprintf(fid,'integer(ip), parameter :: MC_%s_NUM_CASES = %d\n',elem_type,mc_ncases);
fprintf(fid,'integer(ip), parameter :: MC_%s_MAX_NUM_SUBCELLS = %d\n',elem_type,mc_max_sub_cells);
fprintf(fid,'integer(ip), parameter :: MC_%s_NUM_NODES_PER_SUBCELL = %d\n',elem_type,mc_num_nodes_per_subcell);

%fprintf(fid,'\n');

%fprintf(fid,'integer(ip), parameter :: MC_%s_NUM_SUBCELLS_PER_CASE(MC_%s_NUM_CASES) = &\n',elem_type,elem_type);
if mc_ncases < 50
    endaux = '';
else
    endaux = '&\n';
end
fprintf(fid,['integer(ip), parameter :: MC_%s_NUM_SUBCELLS_PER_CASE(%d) = ' endaux],elem_type,mc_ncases);

switch elem_type
    case 'HEX8'        
        width = 16;        
    otherwise
        width = mc_ncases;
end

write_long_vector(fid,mc_num_sub_cells_per_case,width,0);
fprintf(fid,'\n');

fprintf(fid,'\n');

fprintf(fid,['integer(ip), parameter :: MC_%s_NUM_CUT_EDGES_PER_CASE(%d) = ' endaux],elem_type,mc_ncases);
write_long_vector(fid,mc_num_cut_edges_per_case,width,0);
fprintf(fid,'\n');

%fprintf(fid,'\n');

switch elem_type
    case 'HEX8'        
        width = mc_max_sub_cells*mc_num_nodes_per_subcell/2;        
    otherwise
        width = mc_max_sub_cells*mc_num_nodes_per_subcell;
end

%fprintf(fid,'integer(ip), parameter :: MC_%s_SUBCELLS_PER_CASE(MC_%s_NUM_NODES_SUBCELL,MC_%s_MAX_NUM_SUBCELLS,MC_%s_NUM_CASES) = reshape( &\n',elem_type,elem_type,elem_type,elem_type);
fprintf(fid,'integer(ip), parameter :: MC_%s_SUBCELL_NODE_IDS_PER_CASE(%d,%d,%d) = &\nreshape( ',elem_type,mc_num_nodes_per_subcell,mc_max_sub_cells,mc_ncases);
write_long_vector(fid,permute(mc_subcells_per_case,[3 2 1]),width,10);
fprintf(fid,' , [%d,%d,%d] )\n',mc_num_nodes_per_subcell,mc_max_sub_cells,mc_ncases);

%fprintf(fid,'\n');

%fprintf(fid,'integer(ip), parameter :: MC_%s_INOUT_SUBCELLS_PER_CASE(MC_%s_MAX_NUM_SUBCELLS,MC_%s_NUM_CASES) = reshape( &\n',elem_type,elem_type,elem_type);
fprintf(fid,'integer(ip), parameter :: MC_%s_INOUT_SUBCELLS_PER_CASE(%d,%d) = &\nreshape( ',elem_type,mc_max_sub_cells,mc_ncases);
write_long_vector(fid,mc_inout_subcells_per_case',mc_max_sub_cells,10);
fprintf(fid,' , [%d,%d] )\n',mc_max_sub_cells,mc_ncases);



fc = fclose(fid);


function write_long_vector(fid,data,width,skip)



cell_bunch = AUX_compute_cell_bunch(numel(data),width);

if isempty(cell_bunch{end})
    cell_bunch = cell_bunch(1:end-1);
end

for ic = 1:length(cell_bunch)
    
    if ic > 1
        fprintf(fid,['%' num2str(skip) 's'],' ');
    else
        fprintf(fid,'[');
    end
        
    for jc = 1:length(cell_bunch{ic})
        if ic == length(cell_bunch) &&  jc == length(cell_bunch{ic})
            str_aux = '%2d';
        else
            str_aux = '%2d,';
        end
        fprintf(fid,str_aux,data(cell_bunch{ic}(jc)));
    end
    
    if ic == length(cell_bunch)
        fprintf(fid,' ]');
    else
        fprintf(fid,' &\n');
    end
end

function cell_bunch = AUX_compute_cell_bunch(nbel,nbel_bunch)

n = floor(nbel/nbel_bunch);
cell_bunch = cell(n+1,1);
aux_1 = 1:nbel_bunch;
for i = 1:n+1
    if i <=n
        cell_bunch{i} =  aux_1+nbel_bunch*(i-1);
    else
        cell_bunch{i} = nbel_bunch*n+1:nbel;
    end
end

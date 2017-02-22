clear
close all
clc

Xe     = [-1 -1; 1 -1; -1 1; 1 1]; % fempar numbering
Eedges = [1 2; 3 4; 1 3; 2 4]; % fempar numbering
Pe  = [0 0 0 0]';

mc_ncases = 2^4;
mc_max_sub_cells = 0;
mc_num_sub_cells_per_case = zeros(mc_ncases,1);
mc_subcells_per_case_aux = cell(mc_ncases,1);
mc_inout_subcells_per_case_aux = cell(mc_ncases,1);
mc_num_nodes_per_subcell = 3;
mc_num_cut_edges_per_case = zeros(mc_ncases,1);


node2bit = [1 2 4 8 16 32 64 128];

Pn = [-1 1];
for n1 = 1:2
    Pe(1)=Pn(n1);
    for n2 = 1:2
        Pe(2)=Pn(n2);
        for n3 = 1:2
            Pe(3)=Pn(n3);
            for n4 = 1:2
                Pe(4)=Pn(n4);
                %disp(icase)
                %disp(Pe')
                
                icase = 0;
                for i=1:length(Pe)
                    if Pe(i) <0
                    icase = bitor(icase,node2bit(i));
                    end
                end
                icase = icase + 1;
                
                [Xtris,Ttris,Ptris] = subtriangulate_element(Xe,Pe,Eedges);
                num_sub_cells = size(Ttris,1);
                
                mc_max_sub_cells = max([mc_max_sub_cells num_sub_cells]);
                mc_num_sub_cells_per_case(icase) = num_sub_cells;
                mc_subcells_per_case_aux{icase} = Ttris;
                mc_inout_subcells_per_case_aux{icase} = Ptris;
                
                if size(Xtris,1)-size(Xe,1) < 0
                    mc_num_cut_edges_per_case(icase) = 0;
                else
                    mc_num_cut_edges_per_case(icase) = size(Xtris,1)-size(Xe,1);
                end
                
                
            end
        end
    end
end


% Transform cells to arrays, now that we know the maximum number of
% subcells
mc_subcells_per_case = zeros(mc_ncases,mc_max_sub_cells,mc_num_nodes_per_subcell);
mc_inout_subcells_per_case = zeros(mc_ncases,mc_max_sub_cells);
for icase = 1:mc_ncases
    N = mc_num_sub_cells_per_case(icase);
    mc_subcells_per_case(icase,1:N,:) = mc_subcells_per_case_aux{icase};
    mc_inout_subcells_per_case(icase,1:N) = mc_inout_subcells_per_case_aux{icase};
end

mc_max_num_cut_edges = max(mc_num_cut_edges_per_case);


elem_type = 'QUA4';
file_name = '../mc_tables_qua4.i90';
write_i90_file(file_name, mc_ncases, mc_max_sub_cells, mc_max_num_cut_edges, ...
mc_num_sub_cells_per_case, mc_subcells_per_case, mc_inout_subcells_per_case,...
mc_num_nodes_per_subcell,mc_num_cut_edges_per_case,elem_type);

disp('Table for QUA4 done!');




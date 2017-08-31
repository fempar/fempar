clear
close all
clc

%addpath('~/Workspace/code/matlab/bddc/libFEM/')

Xe = [...
    -1.0  -1.0  -1.0
    1.0  -1.0  -1.0
    -1.0   1.0  -1.0
    1.0   1.0  -1.0
    -1.0  -1.0   1.0
    1.0  -1.0   1.0
    -1.0   1.0   1.0
    1.0   1.0   1.0 ]; % Fempar numeration

Eedges = [...
    1   2
    3   4
    5   6
    7   8
    1   3
    2   4
    5   7
    6   8
    1   5
    2   6
    3   7
    4   8]; % Fempar numeration

%Efaces = [1 3 2; 1 2 4; 1 4 3; 2 3 4]; % Matlab Numeration

Efaces =  [1 2 3; 1 4 2; 1 3 4; 2 4 3]; % FEMPAR umeration We have to orinte the faces like the first face in fempar (i.e. pointing inwards)


node2bit = [1 2 4 8 16 32 64 128];

Pe  = zeros(8,1);

mc_ncases = 2^8;
mc_max_sub_cells = 0;
mc_max_sub_faces = 0;
mc_num_sub_cells_per_case = zeros(mc_ncases,1);
mc_num_sub_faces_per_case = zeros(mc_ncases,1);
mc_subcells_per_case_aux = cell(mc_ncases,1);
mc_subfacets_per_case_aux = cell(mc_ncases,1);
mc_inout_subcells_per_case_aux = cell(mc_ncases,1);
mc_num_nodes_per_subcell = 4;
mc_num_nodes_per_subfacet = 3;
mc_num_cut_edges_per_case = zeros(mc_ncases,1);

icase = 1;

Pn = [-1 1];
for n1 = 1:2
    Pe(1)=Pn(n1);
    for n2 = 1:2
        Pe(2)=Pn(n2);
        for n3 = 1:2
            Pe(3)=Pn(n3);
            for n4 = 1:2
                Pe(4)=Pn(n4);
                for n5 = 1:2
                    Pe(5)=Pn(n5);
                    for n6 = 1:2
                        Pe(6)=Pn(n6);
                        for n7 = 1:2
                            Pe(7)=Pn(n7);
                            for n8 = 1:2
                                Pe(8)=Pn(n8);
                                
                                icase = 0;
                                for i=1:length(Pe)
                                    if Pe(i) <0
                                        icase = bitor(icase,node2bit(i));
                                    end
                                end
                                icase = icase + 1;
                                
                                [Xtris,Ttris,Ptris] = subtriangulate_element(Xe,Pe,Eedges);
                                num_sub_cells = size(Ttris,1);
                                
                                
                                
                                
                                
                                % Get the interface faces --------------------------(begin)
                                Ttris_faces = [];
                                
                                % Find neigbours of sub-elements
                                [neig, nface] = FEM_find_face_neighbours_3d(Ttris,Efaces);
                                
                                % Find the interior sub-elements
                                intsubelems = find(Ptris<0);
                                
                                % Take the faces, whose neigbour is outside
                                for iaux = 1:length(intsubelems)
                                    isubelem = intsubelems(iaux);
                                    
                                    % Find faces with neighbours
                                    faces = find( neig(isubelem,:)~=0 );
                                    
                                    % Take the faces, whose neigbour is outside
                                    faces = faces(Ptris(neig(isubelem,faces))>0);
                                    
                                    % Store them
                                    for ifac = 1:length(faces)
                                        Ttris_faces = [Ttris_faces; Ttris(isubelem,Efaces(faces(ifac),:))  ]; %#ok<AGROW>
                                    end
                                    
                                end
                                
                                num_sub_faces = size(Ttris_faces,1);
                                
%                                 figure(1)
%                                 clf
%                                 if not(isempty(Xtris))
%                                     FEM_plot_mesh_3d(Xtris,Ttris(Ptris<0,:),'k','FaceColor','r','FaceAlpha',0.5)
%                                     FEM_plot_faces_3d(Xtris,Ttris_faces,1)
%                                 end
%                                 axis equal off
%                                 colorbar
%                                 pause
                                
                                % ---------------------------------------------------(end)
                                
                                
                                
                                
                                
                                
                                
                                
                                mc_max_sub_cells = max([mc_max_sub_cells num_sub_cells]);
                                mc_max_sub_faces = max([mc_max_sub_faces num_sub_faces]);
                                mc_num_sub_cells_per_case(icase) = num_sub_cells;
                                mc_num_sub_faces_per_case(icase) = num_sub_faces;
                                mc_subcells_per_case_aux{icase} = Ttris;
                                mc_inout_subcells_per_case_aux{icase} = Ptris;
                                mc_subfacets_per_case_aux{icase} = Ttris_faces;
                                
                                if size(Xtris,1)-size(Xe,1) < 0
                                    mc_num_cut_edges_per_case(icase) = 0;
                                else
                                    mc_num_cut_edges_per_case(icase) = size(Xtris,1)-size(Xe,1);
                                end
                                
                                
                                
                            end
                        end
                    end
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

% Transform cells to arrays, now that we know the maximum number of
% subfacets
mc_subfacets_per_case = zeros(mc_ncases,mc_max_sub_faces,mc_num_nodes_per_subfacet);
for icase = 1:mc_ncases
    N = mc_num_sub_faces_per_case(icase);
    mc_subfacets_per_case(icase,1:N,:) = mc_subfacets_per_case_aux{icase};
end


elem_type = 'HEX8';
file_name = '../mc_tables_hex8.i90';
write_i90_file(file_name, mc_ncases, mc_max_sub_cells, mc_max_num_cut_edges, ...
    mc_num_sub_cells_per_case, mc_subcells_per_case, mc_inout_subcells_per_case,...
    mc_num_nodes_per_subcell,mc_num_cut_edges_per_case,elem_type,mc_num_nodes_per_subfacet,...
    mc_max_sub_faces,mc_num_sub_faces_per_case,mc_subfacets_per_case);
disp('Table for HEX8 done!')


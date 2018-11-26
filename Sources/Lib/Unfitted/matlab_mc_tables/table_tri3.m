clear
close all
clc

% addpath('~/Workspace/code/matlab/bddc/libFEM/')

Xe     = [0 0; 1 0; 0 1]; % fempar numbering
Eedges = [1 2; 1 3; 2 3]; % fempar numbering
Pe  = [0 0 0 ]';
Efaces = [1 2; 2 3; 3 1]; % Matlab numbering

mc_ncases = 2^3;
mc_max_sub_cells = 0;
mc_max_sub_faces = 0;
mc_max_sub_vefs = 0;

mc_num_vefs = size(Eedges,1);
mc_num_facets = mc_num_vefs;
mc_num_sub_cells_per_case = zeros(mc_ncases,1);
mc_num_sub_faces_per_case = zeros(mc_ncases,1);
mc_num_sub_vefs_per_case = zeros(mc_ncases,mc_num_vefs);
mc_facet_type_per_case_and_facet = zeros(mc_ncases,mc_num_vefs);
mc_subcells_per_case_aux = cell(mc_ncases,1);
mc_subfacets_per_case_aux = cell(mc_ncases,1);
mc_subvefs_per_case_aux = cell(mc_ncases,1);
mc_subvefs_inout_per_case_aux = cell(mc_ncases,1);
mc_inout_subcells_per_case_aux = cell(mc_ncases,1);
mc_num_nodes_per_subcell = 3;
mc_num_nodes_per_subfacet = 2;
mc_num_cut_edges_per_case = zeros(mc_ncases,1);


node2bit = [1 2 4 8 16 32 64 128];

Pn = [-1 1];
for n1 = 1:2
    Pe(1)=Pn(n1);
    for n2 = 1:2
        Pe(2)=Pn(n2);
        for n3 = 1:2
            Pe(3)=Pn(n3);
                
                
                icase = 0;
                for i=1:length(Pe)
                    if Pe(i) <0
                        icase = bitor(icase,node2bit(i));
                    end
                end
                icase = icase + 1;
                
                %                 disp('========================')
                %                 disp(icase)
                %                 disp(Pe')
                
                [Xtris,Ttris,Ptris] = subtriangulate_element(Xe,Pe,Eedges);
                num_sub_cells = size(Ttris,1);
                
                % Get the interface faces --------------------------(begin)
                Ttris_faces = [];
                
                % Find neigbours of sub-elements
                [neig, nface] = FEM_find_neighbours_2d(Ttris,Efaces);
                
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
                
                [Tf, Tf_neigs] = FEM_find_mesh_faces(Ttris);
                
                Ttris_faces_bound = cell(size(Eedges,1),1);
                Ttris_faces_bound_inout = cell(size(Eedges,1),1);
                
                % Find the faces touching each face
                for iface = 1:size(Eedges,1)
                    
                    Piface = Pe(Eedges(iface,:));
                    
                    if ( Piface(1)*Piface(2) < 0 )
                        mc_facet_type_per_case_and_facet(icase,iface) = 0;
                    elseif  Piface(1) > 0
                        mc_facet_type_per_case_and_facet(icase,iface) = 1;
                    else
                        mc_facet_type_per_case_and_facet(icase,iface) = -1;
                    end
                    
                    Ttris_faces_i = [];
                    Ttris_faces_inout_i = [];
                    
                    num_sub_vefs = 0;
                    
                    if mc_facet_type_per_case_and_facet(icase,iface) == 0
                        
                        v1 = [0 0];
                        
                        v = Xe(Eedges(iface,2),:) - Xe(Eedges(iface,1),:);
                        
                        v1(1) =  v(2);
                        v1(2) = -v(1);
                        
                        
                        
                        for jface = 1:size(Tf,1)
                            
                            Xf = Xtris(Tf(jface,:),:);
                            
                            two_nodes_on = 1;
                            for inod = 1:size(Xf,1)
                                v2 = Xf(inod,:) - Xe(Eedges(iface,1),:);
                                if ( abs(dot(v1,v2)) > 1e-12 )
                                    two_nodes_on = 0;
                                end
                            end
                            
                            if (two_nodes_on == 1)
                                
                                inout = Ptris(Tf_neigs(jface,1));
                                
                                Ttris_faces_inout_i = [Ttris_faces_inout_i; inout];
                                Ttris_faces_i = [Ttris_faces_i; Tf(jface,:)];
                                
                            end
                            
                            num_sub_vefs = size(Ttris_faces_i,1);
                            mc_max_sub_vefs = max([mc_max_sub_vefs num_sub_vefs]);
                            
                        end
                        
                    end
                    
                    mc_num_sub_vefs_per_case(icase,iface) = num_sub_vefs;
                    
                    %                     figure(1)
                    %                     clf
                    %                     FEM_plot_element_wise_constant_field_2d(Ptris,Xtris,Ttris,'EdgeColor','k')
                    %                     FEM_plot_mesh_2d(Xtris,Ttris,'k',1,1)
                    %                     FEM_plot_Neumann_faces_2d(Xtris,Ttris_faces_i(Ttris_faces_inout_i>0,:),'r')
                    %                     FEM_plot_Neumann_faces_2d(Xtris,Ttris_faces_i(Ttris_faces_inout_i<0,:),'b')
                    %                     axis equal off
                    %                     colorbar
                    %                     pause
                    
                    Ttris_faces_bound{iface} = Ttris_faces_i;
                    Ttris_faces_bound_inout{iface} = Ttris_faces_inout_i;
                end
                
                
                
%                 
%                                 figure(1)
%                                 clf
%                                 FEM_plot_element_wise_constant_field_2d(Ptris,Xtris,Ttris,'EdgeColor','k')
%                                 FEM_plot_Neumann_faces_2d(Xtris,Ttris_faces)
%                                 axis equal off
%                                 colorbar
%                                 pause(1)
                
                % ---------------------------------------------------(end)
                
                mc_max_sub_cells = max([mc_max_sub_cells num_sub_cells]);
                mc_max_sub_faces = max([mc_max_sub_faces num_sub_faces]);
                mc_num_sub_cells_per_case(icase) = num_sub_cells;
                mc_num_sub_faces_per_case(icase) = num_sub_faces;
                mc_subcells_per_case_aux{icase} = Ttris;
                mc_inout_subcells_per_case_aux{icase} = Ptris;
                mc_subfacets_per_case_aux{icase} = Ttris_faces;
                mc_subvefs_per_case_aux{icase} = Ttris_faces_bound;
                mc_subvefs_inout_per_case_aux{icase} = Ttris_faces_bound_inout;
                
                if size(Xtris,1)-size(Xe,1) < 0
                    mc_num_cut_edges_per_case(icase) = 0;
                else
                    mc_num_cut_edges_per_case(icase) = size(Xtris,1)-size(Xe,1);
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


% Transform cells to arrays, now that we know the maximum number of
% subfacets

mc_subvefs_per_case = zeros(mc_ncases,mc_num_vefs,mc_max_sub_vefs,mc_num_nodes_per_subfacet);
mc_subvefs_inout_per_case = zeros(mc_ncases,mc_num_vefs,mc_max_sub_vefs);
for icase = 1:mc_ncases
    for iface = 1:mc_num_vefs
        N = mc_num_sub_vefs_per_case(icase,iface);
        mc_subvefs_per_case(icase,iface,1:N,:) = mc_subvefs_per_case_aux{icase}{iface};
        mc_subvefs_inout_per_case(icase,iface,1:N) = mc_subvefs_inout_per_case_aux{icase}{iface};
    end
end

elem_type = 'TRI3';
file_name = '../mc_tables_tri3.i90';
write_i90_file(file_name, mc_ncases, mc_max_sub_cells, mc_max_num_cut_edges, ...
    mc_num_sub_cells_per_case, mc_subcells_per_case, mc_inout_subcells_per_case,...
    mc_num_nodes_per_subcell,mc_num_cut_edges_per_case,elem_type,mc_num_nodes_per_subfacet,...
    mc_max_sub_faces,mc_num_sub_faces_per_case,mc_subfacets_per_case,...
    mc_max_sub_vefs,mc_num_sub_vefs_per_case,mc_subvefs_per_case,mc_subvefs_inout_per_case,mc_num_facets,mc_facet_type_per_case_and_facet);

disp('Table for TRI3 done!');




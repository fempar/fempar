clear
close all
clc


Xe = [...
    0.0  0.0  0.0 
    1.0  0.0  0.0
    0.0  1.0  0.0
    0.0  0.0  1.0 ]; % Fempar numeration

Eedges = [...
    1   2
    1   3
    2   3
    1   4
    2   4
    3   4 ]; % Fempar numeration

Efaces = [...
    1 2 3
    1 2 4
    1 3 4
    2 3 4 ]; % Fempar numeration

%Efaces = [1 3 2; 1 2 4; 1 4 3; 2 3 4]; % Matlab Numeration

Efaces_subelem =  [1 2 3; 1 4 2; 1 3 4; 2 4 3]; % FEMPAR numeration We have to orient the faces like the first face in fempar (i.e. pointing inwards)


node2bit = [1 2 4 8 16 32 64 128];

Pe  = zeros(4,1);

mc_ncases = 2^4;
mc_max_sub_cells = 0;
mc_max_sub_faces = 0;
mc_max_sub_vefs = 0;

mc_num_facets = size(Efaces,1);
mc_num_sub_cells_per_case = zeros(mc_ncases,1);
mc_num_sub_faces_per_case = zeros(mc_ncases,1);
mc_num_sub_vefs_per_case = zeros(mc_ncases,mc_num_facets);
mc_facet_type_per_case_and_facet = zeros(mc_ncases,mc_num_facets);

mc_subcells_per_case_aux = cell(mc_ncases,1);
mc_subfacets_per_case_aux = cell(mc_ncases,1);
mc_subvefs_per_case_aux = cell(mc_ncases,1);
mc_subvefs_inout_per_case_aux = cell(mc_ncases,1);

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
                [neig, nface] = FEM_find_face_neighbours_3d(Ttris,Efaces_subelem);
                
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
                        Ttris_faces = [Ttris_faces; Ttris(isubelem,Efaces_subelem(faces(ifac),:))  ]; %#ok<AGROW>
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
                
                [Tf, Tf_neigs] = FEM_find_faces_3d(Ttris,'tet4');
                %icase
                if (icase ==-1)
                    
                    for kkface = 1:size(Tf,1)
                        figure(1)
                        clf
                        hold on
                        FEM_plot_tet4_mesh(Xtris,Ttris,'k','FaceColor','none')
                        for inod = 1:4
                            text(Xe(inod,1),Xe(inod,2),Xe(inod,3),int2str(inod))
                        end
                        Tf(kkface,:)
                        if not(isempty(find(Tf(kkface,:)==5, 1)))
                            Tf(kkface,:)
                            FEM_plot_faces_3d(Xtris,Tf(kkface,:),'r')
                            FEM_plot_tet4_mesh(Xtris,Ttris(Tf_neigs(kkface,1),:),'k','FaceColor','red','FaceAlpha',0.1)
                            if Tf_neigs(kkface,2) ~=0
                                FEM_plot_tet4_mesh(Xtris,Ttris(Tf_neigs(kkface,2),:),'k','FaceColor','blue','FaceAlpha',0.1)
                            end
                            axis equal off
                            pause
                        end
                        
                    end
                    
                end
                
                
                
                
                % Find the faces touching each face (Loop in big faces)
                for iface = 1:size(Efaces,1)
                    
                    is_face_cut = 0;
                    Piface = Pe(Efaces(iface,:));
                    has_neg = not(isempty(find(Piface<0, 1)));
                    has_pos = not(isempty(find(Piface>0, 1)));
                    if (has_neg && has_pos)
                        mc_facet_type_per_case_and_facet(icase,iface) = 0;
                    elseif has_neg
                        mc_facet_type_per_case_and_facet(icase,iface) = -1;
                    else
                        mc_facet_type_per_case_and_facet(icase,iface) = 1;
                    end
                    
                    Ttris_faces_i = [];
                    Ttris_faces_inout_i = [];
                    num_sub_vefs = 0;
                    
                    if (mc_facet_type_per_case_and_facet(icase,iface) == 0)
                        
                        % Compute normal vector big face
                        a1 = Xe(Efaces(iface,2),:) - Xe(Efaces(iface,1),:);
                        a2 = Xe(Efaces(iface,3),:) - Xe(Efaces(iface,1),:);
                        v1 = cross(a1,a2);
                        
                        for jface = 1:size(Tf,1)
                            
                            Xf = Xtris(Tf(jface,:),:);
                            two_three_on = 1;
                            
                            for inod = 1:size(Xf,1)
                                v2 = Xf(inod,:) - Xe(Efaces(iface,1),:);
                                if ( abs(dot(v1,v2)) > 1e-9 && not(norm(v2)<1e-5) )
                                    two_three_on = 0;
                                end
                            end
                            
                            if (two_three_on == 1)
                                inout = Ptris(Tf_neigs(jface,1));
                                Ttris_faces_inout_i = [Ttris_faces_inout_i; inout];
                                Ttris_faces_i = [Ttris_faces_i; Tf(jface,:)];
                            end
                            
                            num_sub_vefs = size(Ttris_faces_i,1);
                            mc_max_sub_vefs = max([mc_max_sub_vefs num_sub_vefs]);
                            
                        end
                        
                    end
                    
                    mc_num_sub_vefs_per_case(icase,iface) = num_sub_vefs;
                    Ttris_faces_bound{iface} = Ttris_faces_i;
                    Ttris_faces_bound_inout{iface} = Ttris_faces_inout_i;
                    
                    
                    if (icase == -1)
                        Piface
                        mc_facet_type_per_case_and_facet(icase,iface)
                        figure(1)
                        clf
                        hold on
                        FEM_plot_tet4_mesh(Xtris,Ttris(Ptris<0,:),'k','FaceColor','red','FaceAlpha',0.1)
                        FEM_plot_tet4_mesh(Xtris,Ttris(Ptris>0,:),'k','FaceColor','blue','FaceAlpha',0.1)
                        FEM_plot_faces_3d(Xtris,Ttris_faces_i(Ttris_faces_inout_i>0,:),1,'b')
                        FEM_plot_faces_3d(Xtris,Ttris_faces_i(Ttris_faces_inout_i<0,:),1,'r')
                        for inod = 1:4
                            text(Xe(inod,1),Xe(inod,2),Xe(inod,3),int2str(inod))
                        end
                        
                        axis equal off
                        colorbar
                        pause
                        
                    end
                    
                    
                    
                end
                
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

mc_subvefs_per_case = zeros(mc_ncases,mc_num_facets,mc_max_sub_vefs,mc_num_nodes_per_subfacet);
mc_subvefs_inout_per_case = zeros(mc_ncases,mc_num_facets,mc_max_sub_vefs);
for icase = 1:mc_ncases
    for iface = 1:mc_num_facets
        N = mc_num_sub_vefs_per_case(icase,iface);
        mc_subvefs_per_case(icase,iface,1:N,:) = mc_subvefs_per_case_aux{icase}{iface};
        mc_subvefs_inout_per_case(icase,iface,1:N) = mc_subvefs_inout_per_case_aux{icase}{iface};
    end
end


elem_type = 'TET4';
file_name = '../mc_tables_tet4.i90';
write_i90_file(file_name, mc_ncases, mc_max_sub_cells, mc_max_num_cut_edges, ...
    mc_num_sub_cells_per_case, mc_subcells_per_case, mc_inout_subcells_per_case,...
    mc_num_nodes_per_subcell,mc_num_cut_edges_per_case,elem_type,mc_num_nodes_per_subfacet,...
    mc_max_sub_faces,mc_num_sub_faces_per_case,mc_subfacets_per_case,...
    mc_max_sub_vefs,mc_num_sub_vefs_per_case,mc_subvefs_per_case,mc_subvefs_inout_per_case,mc_num_facets,mc_facet_type_per_case_and_facet);
disp('Table for TET4 done!')


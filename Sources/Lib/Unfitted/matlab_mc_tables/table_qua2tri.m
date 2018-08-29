clear
close all
clc

quad_faces       = [ 1 2; 3 4; 1 3 ; 2 4];
quad_coordinates = [-1 -1; 1 -1; -1 1; 1 1];

tri_faces = [1 2; 1 3; 2 3];
tri_coordinates = [0 0; 1 0; 0 1]



tri_cells       = delaunay(quad_coordinates);
num_tri_cells 	= size(tri_cells, 1);
nodes_per_tri_cell = size(tri_cells,2);


%% Reorder according to the FEMPAR numeration of triangles 
%%  
%%  (0,1)
%%    3
%%    |`.
%%    |  `. (3) 
%% (2)|    `.
%%    |      `.  
%%    |________`.
%%    1   (1)   2
%%  (0,0)     (1,0)

origin_vertex_num = 1;

for icell=1:num_tri_cells
	for inode=1:nodes_per_tri_cell
		aux = tri_cells;
		aux(icell,:) = [];
		% Check shared is the node is shared with other cells, i.e., origin vertex
		if(isempty(find(aux==tri_cells(icell,inode))) && inode~=origin_vertex_num)
			tri_cells(icell,:) = circshift(tri_cells(icell,:),origin_vertex_num-inode)
		end
	end
end 




num_nodes_per_quad_face = size(quad_faces,2);
num_nodes_per_trian_face = size(tri_faces,2);


num_quad_faces = size(quad_faces, 1);


for iface=1:num_quad_faces
	id_matching_nodes = ismember(tri_cells,quad_faces(iface,:));
	% max(...): limits to one entry when one quad face have more than one triangle face
    face_in_quad_to_trian(iface) = max( find( sum(id_matching_nodes,2) == num_nodes_per_trian_face ) );

    % Get the triangle nodes in local FEMPAR numbering
    local_trian_face = find( id_matching_nodes(face_in_quad_to_trian(iface),:)) ; 
    % Match face numbering with list in 'tri_faces' to get local face id
    face_in_quad_to_face_in_trian(iface) = find(sum(ismember(tri_faces, local_trian_face),2) == num_nodes_per_trian_face);
end


file_name = '../mc_tables_quad2trian.i90';

write_conversion_i90_file(file_name, num_quad_faces, face_in_quad_to_trian, face_in_quad_to_face_in_trian )

disp('Table for conversion QUA4 to TRI3 done!');
% face_in_quad_to_trian
% face_in_quad_to face_in_trian



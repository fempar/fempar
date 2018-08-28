clear
close all
clc

quad_faces       = [ 1 2; 3 4; 1 3 ; 2 4];
quad_coordinates = [-1 -1; 1 -1; -1 1; 1 1];

tri_faces = [1 2; 2 3; 3 1];



tri_cells       = delaunay(quad_coordinates);
num_tri_cells 	= size(tri_cells, 1);
nodes_per_tri_cell = size(tri_cells,2);


%% Reorder according to the FEMPAR numeration of triangles (GUESS)
%%
%% (2)__3__(3)
%%    \    |
%%     \   |
%%     1\  |2
%%       \ |
%%        \|
%%        (1)
%%

vertex_num = 3;

for icell=1:num_tri_cells
	for inode=1:nodes_per_tri_cell
		aux = tri_cells;
		aux(icell,:) = [];
		if(isempty(find(aux==tri_cells(icell,inode))) && inode~=vertex_num)
			tri_cells(icell,:) = circshift(tri_cells(icell,:),vertex_num-inode)
		end
	end
end 






num_nodes_per_quad_face = 2;
num_nodes_per_trian_face = 2;


num_quad_faces = size(quad_faces, 1);


for iface=1:num_quad_faces
	for icell=1:num_tri_cells
        idnodes = -1*ones(1,num_nodes_per_quad_face);
		for inode=1:num_nodes_per_quad_face
			id = find(tri_cells(icell,:)==quad_faces(iface,inode));
			if ( ~isempty(id) )
				idnodes(inode) = id;	
			end
		end
		if( min(idnodes)>-1)
			face_in_quad_to_trian(iface) = icell;
			for inode=1:length(idnodes)
				imax = inode+num_nodes_per_trian_face-1; 
				irange = inode:imax;
				if ( imax > length(idnodes) )
					irange = [inode:length(idnodes),1:imax-length(idnodes)];
				end
				id_aux=idnodes(irange);

				if( isempty(id_aux(diff(id_aux)~=1)))
					face_in_quad_to_face_in_trian(iface) = min(id_aux);
				elseif ( length(id_aux(diff(id_aux)~=1)) == 1 && isempty(id_aux(diff(id_aux)==-1)) )
					face_in_quad_to_face_in_trian(iface) = min(id_aux(find((diff(id_aux)~=1)==1)+1:end));
				end
			end
		end

	end
end

file_name = '../mc_tables_quad2trian.i90';

write_conversion_i90_file(file_name, num_quad_faces, face_in_quad_to_trian, face_in_quad_to_face_in_trian )

disp('Table for conversion QUA4 to TRI3 done!');
% face_in_quad_to_trian
% face_in_quad_to face_in_trian



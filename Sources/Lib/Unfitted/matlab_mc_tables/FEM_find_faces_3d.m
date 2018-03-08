function [faces_nodes,faces_neig,faces_neig_faces] = FEM_find_faces_3d(T,elem_type)

% Var              Size                        Description
% ================ =========================== ==================================
% T                [nbel,nben]
% elem_type        str
% faces_nodes      [nbfac,nbfn]
% faces_neig       [nbfac,2]
% faces_neig_faces [nbfac,2]

switch elem_type
    case 'hex8'
        Faces = [1 4 3 2; 1 2 6 5; 2 3 7 6; 3 4 8 7;  4 1 5 8 ; 5 6 7 8];
    case 'tet4'
         Faces = [1 3 2; 1 2 4; 2 3 4; 1 4 3];
    otherwise
        disp('ERROR: Element not implemented')
        error('ERROR: Element not implemented')
end

nbel = size(T,1);
nbef = size(Faces,1);
nbfn = size(Faces,2);

% Allocate memory (too much)
faces_nodes      = zeros(nbel*nbef,nbfn);
faces_neig       = zeros(nbel*nbef,2);
faces_neig_faces = zeros(nbel*nbef,2);

% Auxiliary
elems_face_taken = false(nbel,nbef);
ncount = 1;
for iele = 1:nbel
    Te = T(iele,:);
    for ifac = 1:nbef
        if not(elems_face_taken(iele,ifac))
            Tf = Te(Faces(ifac,:));
            jele = FEM_find_neigbour_to_face_3d(T,Tf,iele);
            jfac = FEM_find_neigbour_face_3d(T,Tf,jele,Faces);
            faces_nodes(ncount,:)      = Tf;
            faces_neig(ncount,:)       = [iele jele];
            faces_neig_faces(ncount,:) = [ifac jfac];
            ncount = ncount + 1;
            elems_face_taken(iele,ifac) = true(1);
            if jele ~= 0
                elems_face_taken(jele,jfac) = true(1);
            end
        end
    end
end

% Take rid of empty memory
pos_used = 1:ncount-1;
faces_nodes      = faces_nodes(pos_used,:);
faces_neig       = faces_neig(pos_used,:);
faces_neig_faces = faces_neig_faces(pos_used,:);










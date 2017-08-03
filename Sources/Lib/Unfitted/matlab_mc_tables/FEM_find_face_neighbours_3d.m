function [neig,  neig_face, neig_rota] = FEM_find_face_neighbours_3d(T,Faces)

% Should work for hexs and tets

% Sizes
nbel = size(T,1);
nbef = size(Faces,1);

% Allocate
neig      = zeros(nbel,nbef);
neig_face = zeros(nbel,nbef);
neig_rota = zeros(nbel,nbef);


% Loop in elements
for iel = 1:nbel
    
    Te = T(iel,:);
    
    % Loop in faces
    for ifac = 1:nbef
        
        % Extract nodes of this face
        Tf = Te(Faces(ifac,:));     
        
        % Find Neigbour to that face
        neig(iel,ifac) = Find_neigbour(T,Tf,iel);
        
        % Find Neigbours face and rotation index
        [neig_face(iel,ifac), neig_rota(iel,ifac)] =...
            Find_neigbour_face(T,Tf,neig(iel,ifac),Faces); 
        
    end
end

function neig = Find_neigbour(T,Tf,iel)

nbfn = length(Tf);

% Find the first node agains the rest of the mesh
[neig_cand,~] = find(T == Tf(1));
neig_cand     = setdiff(neig_cand,iel);
T_cand        = T(neig_cand,:);

% Find the other nodes
for inod = 2:nbfn
    [neig_cand_2,~] = find(T_cand == Tf(inod));
    neig_cand = neig_cand(neig_cand_2);
    T_cand    = T(neig_cand,:);
end

% Store neigbour
if not(isempty(T_cand))
    neig = neig_cand;    
else
    neig = 0;    
end

function [neig_face, neig_rota] = Find_neigbour_face(T,Tf,neig,Faces)

if neig == 0
    neig_face = 0;
    neig_rota = 0;
    return
end

% Information about neigbours faces
T_neig = T(neig,:);
T_faces_neig = T_neig(Faces);

% Loop in face nodes
nbfn = length(Tf);
nbef = size(Faces,1);
faces_cand   = 1:nbef;
T_faces_cand = T_faces_neig(faces_cand,:);
for ifn = 1:nbfn    
    [faces_cand_2,~] = find(T_faces_cand == Tf(ifn));    
    faces_cand = faces_cand(faces_cand_2);
    T_faces_cand = T_faces_neig(faces_cand,:);
end

% Set neighbours face
neig_face = faces_cand;

% Find rotation status
Tf_neig = T_faces_neig(neig_face,:);
neig_rota = find(Tf_neig==Tf(1));



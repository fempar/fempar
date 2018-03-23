function [neig_face, neig_rota] = FEM_find_neigbour_face_3d(T,Tf,neig,Faces)

% Var       Size                        Description
% ========= =========================== ==================================
% T         [nbel,nben]
% Tf        [1,nbfn]
% neig      [1,1]
% Faces     [nbef,nbfn]


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
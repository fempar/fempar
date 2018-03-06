function neig = FEM_find_neigbour_to_face_3d(T,Tf,iel)

% Var       Size                        Description
% ========= =========================== ==================================
% T         [nbel,nben]
% Tf        [1,nbfn]
% iel       [1,1]


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

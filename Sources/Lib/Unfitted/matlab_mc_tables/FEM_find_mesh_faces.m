function [Tf,neigs, neigs_f, node_a] = FEM_find_mesh_faces(T)

%
% In
%
% T    (nbel,nben) : Mesh connectivities
%
% Out
%
% Tf      (nbfa,2) : Face connectivities
% neigs   (nbfa,2) : Right and left elements of the faces (not left and right!)
% neigs_f (nbfa,2) : Local face numeration of neigs
% node_a (nbfa,2) : local numerotation of the first node in the oriented face
%


[nbel, nben] = size(T);
nbef = nben;

Efaces = zeros(nbef,2);
Efaces(:,1) = 1:nben;
Efaces(:,2) = [2:nben 1];

Tf      = zeros(0,2);
neigs   = zeros(0,2);
neigs_f = zeros(0,2);
node_a = zeros(0,2);

markE = zeros(nbel,nben);

for iele = 1:nbel
    for iface = 1:nbef
        if not(markE(iele,iface))
            
            [jele, jface] = Find_neighbour(iele,iface,T);
            
            Tf      = [Tf;      T(iele,Efaces(iface,:)) ];
            neigs   = [neigs;   [iele  jele]    ];
            neigs_f = [neigs_f; [iface jface]  ];
            
            markE(iele,iface) = 1;
            if jele ~= 0
                markE(jele,jface) = 1;
            end
            
            
            
            T_aux=[T(iele,:) T(iele,1)];
            
            if iface~=0 && jface~=0
                
                if T_aux([iface] )==Tf(end,1)
                    pos1=1; pos2=2;
                else
                    pos1=2; pos2=1;
                end
                node_a = [node_a; [pos1 pos2] ];
                
                
            elseif iface==0
                if T_aux([ jface] )==Tf(end,1)
                    pos2=2;
                else pos2=1;
                end
                node_a = [node_a; [0 pos2] ];
                
                
            elseif jface==0
                if T_aux([ iface] )==Tf(end,1)
                    pos1=1;
                else pos1=2;
                end
                node_a = [node_a; [pos1 0]];
                
            end
        end
    end
end



function [jele, jface] = Find_neighbour(iele,iface,T)

% In
%
%  iele  (1,1):  Number of element
%  iface (1,1):  Local number of face in iele
%
% Out
%
%  jele  (1,1): Element that touches element iele by the local face iface
%               of iele
%  jface (1,1): Local face of jele that touches iele
%
%  It works just for quadrilaterals
%

nben = size(T,2);
nbef  = nben;

Efaces = zeros(nbef,2);
Efaces(:,1) = 1:nben;
Efaces(:,2) = [2:nben 1];

T_faces = zeros(nbef,2);


nodes_iface = T(iele,Efaces(iface,:));


% Find neighbour

[elems,aux] = find(T==nodes_iface(1));

elems = elems(elems~=iele);
T_elems = T(elems,:);

if (~isempty(elems))
    
    [aux,aux2] = find(T_elems==nodes_iface(2));
    elems = elems(aux);
    
end

if isempty(elems)
    jele = 0;
    jface = 0;
    
else
    jele = elems(1);
    
    % Find local face of the neighbour
    T_faces(:) = T(jele,Efaces);
    [faces aux] = find(T_faces == nodes_iface(1));
    T_faces = T_faces(faces,:);
    [aux aux2] = find(T_faces == nodes_iface(2));
    jface = faces(aux);
    
end





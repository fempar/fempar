function [neig, nface] = FEM_find_neighbours_2d(T,Tfaces)


%
% In
%
% T     (nbel,nben) : Element connectivities
%
% Out
%
% neig  (nbel,nbef) : Neighbours at each face of the elements
% nface (nbel,nbef) : Local face of the neighbour that touches the element
%

[nbel, nben] = size(T);
nbef = nben;


neig = zeros(nbel,nbef);
nface = zeros(nbel,nbef);


for iele = 1:nbel    
    
    for iface = 1:nbef
        
         [jele, jface] = Find_neighbour(iele,iface,T,Tfaces);
        
        neig(iele,iface)  = jele;
        nface(iele,iface) = jface;
        
    end
    
end


function [jele, jface] = Find_neighbour(iele,iface,T,Efaces)

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

T_faces = zeros(nbef,2);


nodes_iface = T(iele,Efaces(iface,:));


% Find neighbour

[elems,~] = find(T==nodes_iface(1));

elems = elems(elems~=iele);
T_elems = T(elems,:);

if (~isempty(elems))

    [aux,~] = find(T_elems==nodes_iface(2));
    elems = elems(aux);

end

if isempty(elems)
    jele = 0;
    jface = 0;

else
    jele = elems(1);

    % Find local face of the neighbour
    T_faces(:) = T(jele,Efaces);
    [faces, ~] = find(T_faces == nodes_iface(1));
    T_faces = T_faces(faces,:);
    [aux, ~] = find(T_faces == nodes_iface(2));
    jface = faces(aux);

end







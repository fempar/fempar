function [Xtris,Ttris,Ptris] = subtriangulate_element(Xe,Pe,Eedges)

global SUBELEMS_INSIDE_BY_DEFAULT

nben = size(Xe,1);
% Chech that we have a cuted element
if length(find(Pe>=0)) == nben || length(find(Pe<=0)) == nben    
    Xtris = [];
    Ttris = [];
    Ptris = [];
    return
end

% We slightly perturbate the zero values to ensure unique delaunay
Pmax = max(abs(Pe));
tol= 1e-12;
Pe(Pe==0) = tol*Pmax;

% Compute the intersection points
nbsd = size(Xe,2);
X0 = zeros(0,nbsd);
X0aux = zeros(0,nbsd);
for ifac = 1:size(Eedges,1)
    nods = Eedges(ifac,:);
    Pnods = Pe(nods);
    if Pnods(1)*Pnods(2) < 0
        Xm = Xe(nods(1),:) + (abs(Pnods(1))./(sum(abs(Pnods),1)))*( Xe(nods(2),:) - Xe(nods(1),:) );
        X0 = [X0; Xm]; %#ok<AGROW>
        X0aux = [X0aux; 0.5*(Xe(nods(1),:) + Xe(nods(2),:)) ];%#ok<AGROW>
    end
end

% Compute trinagultion
Xtris = [Xe;X0];
%Ttris = delaunay(Xtris);
Ttris = delaunay([Xe;X0aux]);

Ttris = reorient_elems(Ttris,[Xe;X0]);


if SUBELEMS_INSIDE_BY_DEFAULT
    
    Ptris = -ones(size(Ttris,1),1); % By default elements are inside
    npos = find(Pe>0)';
    for inod = npos % Only outside the elems that have at least a node outside
        [elems,~] = find(Ttris==inod);
        Ptris(elems) = 1;
    end
    
else
    
    Ptris = ones(size(Ttris,1),1); % By default elements are outside
    nneg = find(Pe<0)';
    for inod = nneg % Only inside the elems that have at least a node inside
        [elems,~] = find(Ttris==inod);
        Ptris(elems) = -1;
    end
    
end

function T = reorient_elems(T,X)

nbel = size(T,1);
nben = size(T,2);

if nben == 3
    return
end

for iel =1:nbel
    
    Te = T(iel,:);
    Xe = X(Te,:);
    
    v1 = Xe(2,:) - Xe(1,:);
    v2 = Xe(3,:) - Xe(1,:);
    v3 = Xe(4,:) - Xe(1,:);
    
    if  dot(v3,cross(v1,v2)) < 0
        Te = Te([1 3 2 4]);
        T(iel,:) = Te;
    end
        
end






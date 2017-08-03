function [D] = edgeLengths(F,V)
%{
compute weighted adjacency matrix between for mesh represnted as shared
vertex data structure

Input - 
        F - faces
        V - vertices

D     - weighted adjacency matrix. the weights are edge lengths


%}

nv = size(V,1);
nf = size(F,1);

% edge length

i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);

l = [ ...
    sqrt(sum((V(i2,:)-V(i3,:)).^2,2)) ...
    sqrt(sum((V(i3,:)-V(i1,:)).^2,2)) ...
    sqrt(sum((V(i1,:)-V(i2,:)).^2,2)) ...
    ];

l1 = l(:,1);
l2 = l(:,2);
l3 = l(:,3);

i = [i1 i2 i2 i3 i3 i1 ];
j = [i2 i1 i3 i2 i1 i3 ];

% values corresponding to pairs form (i,j)
v = [l3 l3 l1 l1 l2 l2];

[C,ia,ic] = unique([i(:),j(:)],'rows');

% distance matrix

%D = sparse(i(:),j(:),v(:),nv,nv);
D = sparse(C(:,1),C(:,2),v(ia),nv,nv);

%D = accumarray([i(:),j(:)],v(:),[nv,nv],@min,0,true);

% ensure symmetry
D = (D+D')/2;

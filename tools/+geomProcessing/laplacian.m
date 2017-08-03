function [W,A] = laplacian(V,F)
% matlab implementation 
% based on A. Jacobson et al.

nv = size(V,1);
nf = size(F,1);

% edge length
l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];

i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);

l1 = l(:,1);
l2 = l(:,2);
l3 = l(:,3);

% Heron's formula
% semiperimeter
s  = (l1 + l2 + l3)*0.5;
% triangle-wise area
fA  = sqrt( s.*(s-l1).*(s-l2).*(s-l3));
    
    
% cotangent weight  
cot12 = (l1.^2 + l2.^2 - l3.^2)./fA/4;
cot23 = (l2.^2 + l3.^2 - l1.^2)./fA/4;
cot31 = (l1.^2 + l3.^2 - l2.^2)./fA/4;
  
diag1 = -cot12-cot31; 
diag2 = -cot12-cot23; 
diag3 = -cot31-cot23;


i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
% values corresponding to pairs form (i,j)
v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];

% stiffness matrix

W = sparse(i,j,v,nv,nv);

% ensure symmetry

%W = (W+W')/2;

% compute vertex-wise area elements
tri2ver = sparse(F(:,1),[1:nf]',ones(nf,1),nv,nf) + ...
          sparse(F(:,2),[1:nf]',ones(nf,1),nv,nf) + ...
          sparse(F(:,3),[1:nf]',ones(nf,1),nv,nf);           
tri2ver(tri2ver>0) = 1;

% barycentric mass matrix
A = spdiags(tri2ver*fA/3,0,nv,nv);

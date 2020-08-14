%”¶±‰æÿ’Û
%181215
function [ B,detJ ] = StraB( ielem,s,t,nodeC,elemN,nnode )
B = zeros(3,2*nnode);
[J1,J] = JacoJ(nodeC,elemN,ielem,nnode,s,t);
detJ = det(J);
temp = J^(-1)*J1;
for i = 1:1:nnode
    B(1,2*(i-1)+1) = temp(1,i);
    B(2,2*(i-1)+2) = temp(2,i);
    B(3,2*(i-1)+1) = temp(2,i);B(3,2*(i-1)+2) = temp(1,i);
end
end


% xiezhuoyu
% mechanics_xzy@163.com
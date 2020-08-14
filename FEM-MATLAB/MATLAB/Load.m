%计算节点荷载列阵
%181215
function [ loadP ] = Load( npoin,lodpt )
loadP = zeros(2*npoin,1);%列阵容量
loadP( (lodpt(1:1:size(lodpt,1),1)-1)*2 + 1 ) = lodpt(1:1:size(lodpt,1),2);%x方向载荷
loadP( (lodpt(1:1:size(lodpt,1),1)-1)*2 + 2 ) = lodpt(1:1:size(lodpt,1),3);%y方向载荷
% for i = 1:1:size(lodpt,1)
%     loadP( (lodpt(i,1)-1)*2 + 1 ) = lodpt(i,2);%x方向载荷
%     loadP( (lodpt(i,1)-1)*2 + 2 ) = lodpt(i,3);%y方向载荷
% end
end


% xiezhuoyu
% mechanics_xzy@163.com
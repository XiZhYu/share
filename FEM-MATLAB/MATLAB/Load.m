%����ڵ��������
%181215
function [ loadP ] = Load( npoin,lodpt )
loadP = zeros(2*npoin,1);%��������
loadP( (lodpt(1:1:size(lodpt,1),1)-1)*2 + 1 ) = lodpt(1:1:size(lodpt,1),2);%x�����غ�
loadP( (lodpt(1:1:size(lodpt,1),1)-1)*2 + 2 ) = lodpt(1:1:size(lodpt,1),3);%y�����غ�
% for i = 1:1:size(lodpt,1)
%     loadP( (lodpt(i,1)-1)*2 + 1 ) = lodpt(i,2);%x�����غ�
%     loadP( (lodpt(i,1)-1)*2 + 2 ) = lodpt(i,3);%y�����غ�
% end
end


% xiezhuoyu
% mechanics_xzy@163.com
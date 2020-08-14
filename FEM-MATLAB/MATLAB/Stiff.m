% ��������նȾ���
% 181215
function [ stifK ] = Stiff( npoin,nelem,nvfix,ntype,nnode,nmats,nhamm,noutp,elemM,elemN, ...
    nodeC,fixN,mater,model,iplod,lodpt,posiH,weigH,DD )
%--------------------------------------------------------------------------
%�鼯����նȾ���
rowK = zeros((2*nnode)^2*nelem,1);%������ռ�洢
colK = zeros((2*nnode)^2*nelem,1);
valK = zeros((2*nnode)^2*nelem,1);
tempJS = 0;%��������
tempA = [ 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9];%һ������Ӧx��y��������
tempB = [ 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];%��ż����Ź���һ��
for ielem = 1:1:nelem
    D = DD(:,:,elemM(ielem));%���Ծ���
    stifeK = Stife(ielem,posiH,weigH,D,nodeC,elemN,nnode);%���㵥Ԫ�նȾ���
    for i = 1:1:2*nnode
        for j = 1:1:2*nnode
            tempJS = tempJS+1;
            rowK(tempJS) = ( elemN(ielem,tempA(i))-1 )*2 + tempB(i);%��¼�к�
            colK(tempJS) = ( elemN(ielem,tempA(j))-1 )*2 + tempB(j);%��¼�к�
            valK(tempJS) = stifeK(i,j);%��¼Ԫ��ֵ
        end
    end
end
stifK = sparse(rowK,colK,valK,2*npoin,2*npoin);%ϡ�����洢
%���Ż�������ϡ�裬���ҶԳ�
%--------------------------------------------------------------------------
end


% xiezhuoyu
% mechanics_xzy@163.com
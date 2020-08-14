% 计算整体刚度矩阵
% 181215
function [ stifK ] = Stiff( npoin,nelem,nvfix,ntype,nnode,nmats,nhamm,noutp,elemM,elemN, ...
    nodeC,fixN,mater,model,iplod,lodpt,posiH,weigH,DD )
%--------------------------------------------------------------------------
%组集整体刚度矩阵
rowK = zeros((2*nnode)^2*nelem,1);%给出大空间存储
colK = zeros((2*nnode)^2*nelem,1);
valK = zeros((2*nnode)^2*nelem,1);
tempJS = 0;%计数变量
tempA = [ 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9];%一个结点对应x、y两个方向
tempB = [ 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];%奇偶数编号规则不一致
for ielem = 1:1:nelem
    D = DD(:,:,elemM(ielem));%弹性矩阵
    stifeK = Stife(ielem,posiH,weigH,D,nodeC,elemN,nnode);%计算单元刚度矩阵
    for i = 1:1:2*nnode
        for j = 1:1:2*nnode
            tempJS = tempJS+1;
            rowK(tempJS) = ( elemN(ielem,tempA(i))-1 )*2 + tempB(i);%记录行号
            colK(tempJS) = ( elemN(ielem,tempA(j))-1 )*2 + tempB(j);%记录列号
            valK(tempJS) = stifeK(i,j);%记录元素值
        end
    end
end
stifK = sparse(rowK,colK,valK,2*npoin,2*npoin);%稀疏矩阵存储
%可优化：不仅稀疏，而且对称
%--------------------------------------------------------------------------
end


% xiezhuoyu
% mechanics_xzy@163.com
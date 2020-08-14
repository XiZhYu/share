%计算单元刚度矩阵
%181215
function [ stifeK ] = Stife( ielem,posiH,weigH,D,nodeC,elemN,nnode )
weigh = [weigH(1) weigH(2) weigH(2) weigH(2) weigH(3) weigH(3) weigH(3)];%方便走循环
%weigh = [A1 B3 B3 B3 C3 C3 C3]
posit = [1/3 1/3;posiH(1) posiH(1);posiH(1) posiH(2);posiH(2) posiH(1) ...
    ; posiH(3) posiH(3);posiH(3) posiH(4);posiH(4) posiH(3)];%方便走循环
%posit = [1/3 1/3;a a;a b;b a;c c;c d;d c]
%Hammer积分
stifeK = zeros(2*nnode,2*nnode);
for i = 1:1:7
    [ B,detJ ] = StraB(ielem,posit(i,1),posit(i,2),nodeC,elemN,nnode);
    stifeK = stifeK + weigh(i)*B'*D*B*detJ;
end
end


% xiezhuoyu
% mechanics_xzy@163.com
% �κ���
% 2018/12/19
function [ N,L ] = ShapN( nnode,s,t )
L = [ s t ];%�������
L3 = 1-sum(L);
switch nnode
    case {3}
%         %���㼸��ϵ��
%         x = nodeC(elemN(ielem,:),1);%�����ζ�������
%         y = nodeC(elemN(ielem,:),2);
%         a = [ x(2)*y(3)-x(3)*y(2) x(3)*y(1)-x(1)*y(3) x(1)*y(2)-x(2)*y(1) ];
%         b = [ -y(3)+y(2) -y(1)+y(3) -y(2)+y(1) ];
%         c = [ -x(3)+x(2) -x(1)+x(3) -x(2)+x(1) ];
% ��Ŀǰ�����ֻ��Ҫ֪�����ֵ�ľֲ����꣨��������꣩������Ҫ�����������������
        N = [ L(1) 0 L(2) 0 L3 0 ; ...
              0 L(1) 0 L(2) 0 L3 ];
    case {6}
        N = [ L(1)*(2*L(1)-1) 0 L(2)*(2*L(2)-1) 0 L3*(2*L3-1) 0 4*L(1)*L(2) 0 4*L(2)*L3 0 4*L3*L(1) 0 ; ...
              0 L(1)*(2*L(1)-1) 0 L(2)*(2*L(2)-1) 0 L3*(2*L3-1) 0 4*L(1)*L(2) 0 4*L(2)*L3 0 4*L3*L(1) ];
    otherwise
        disp('ShapN ERROR!')
end
end


% xiezhuoyu
% mechanics_xzy@163.com
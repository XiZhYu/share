%Hammer积分表格数据
%181215
function [posiH,weigH] = IntHammer( nhamm )
posiH = zeros(4,1);%a,b,c,d
weigH = zeros(3,1);%A1,B3,C3
switch nhamm
    case{1}
        weigH(1) = 0.5;
        weigH(2) = 0.0;
        weigH(3) = 0.0;
    case{2}
        posiH(1) = 1.0/6.0;
        posiH(2) = 2.0/3.0;
        weigH(1) = 0.0;
        weigH(2) = 1.0/6.0;
        weigH(3) = 0.0;
    case{3}
        posiH(1) = 0.2;
        posiH(2) = 0.6;
        weigH(1) = -27.0/96.0;
        weigH(2) = 25.0/96.0;
        weigH(3) = 0.0;
    case{5}
        posiH(1) = 0.1012865073;
        posiH(2) = 0.7974269853;
        posiH(3) = 0.4701420641;
        posiH(4) = 0.0597158717;
        weigH(1) = 0.1125000000;
        weigH(2) = 0.0629695903;
        weigH(3) = 0.0661970764;
    otherwise
        disp('IntHammer ERROR!')
end
end


% xiezhuoyu
% mechanics_xzy@163.com
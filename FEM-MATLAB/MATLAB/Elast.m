%���Ծ���
%181215
function [ elasD ] = Elast( mater,ntype )
elasD = zeros(3,3,size(mater,1));
for i = 1:1:size(mater,1)
    switch ntype
        case {1}%ƽ��Ӧ��
            E0 = mater(i,2);%E0 = E;
            Mu0 = mater(i,3);%Mu0 = Mu;
        case {2}%ƽ��Ӧ��
            E0 = mater(i,2)/(1-mater(i,3)^2);%E0 = E/(1-Mu^2);
            Mu0 = mater(i,3)/(1-mater(i,3));%Mu0 = Mu/(1-Mu);
        otherwise
            disp('��ƽ��Ӧ����ƽ��Ӧ������')
    end
    elasD(:,:,i) = E0/(1-Mu0^2) * [1 Mu0 0;Mu0 1 0;0 0 (1-Mu0)/2];
end
end

% xiezhuoyu
% mechanics_xzy@163.com

% ���Դ����������ǰ����
% 18/12/16
function [ loadP,stifK ] = PreLAES( loadPf,stifKf,fixN )
method = 1;
switch method
    %������
    case(1)
        LG = 1E+10;%����
        for i = 1:1:size(fixN,1)
            x = (fixN(i,1)-1)*2 + 1;
            switch fixN(i,2)
                case(10)%x����
                    loadPf( x ) = fixN(i,3)*LG;
                    stifKf( x,x ) = LG;
                case(01)%y���򣬵ȼ�1
                    loadPf( x+1 ) = fixN(i,4)*LG;
                    stifKf( x+1,x+1 ) = LG;
                case(11)%x��y����
                    loadPf( x ) = fixN(i,3)*LG;
                    loadPf( x+1 ) = fixN(i,4)*LG;
                    stifKf( x , x ) = LG;
                    stifKf( x+1 , x+1 ) = LG;
                otherwise
                    disp('DEBUG: LargeNumber!')
            end
        end
        loadP = loadPf;
        stifK = stifKf;
    %���з�
    %�ƻ��˾���ԭ�����ݡ���С
    case(2)
        ID = 1E+30;%��Ϊ��ǣ������������ɾ��
        for i = 1:1:size(fixN,1)
            x = (fixN(i,1)-1)*2 + 1;
            switch fixN(i,2)
                case(10)%x����
                    loadPf( x ) = ID;
                    stifKf( x,x ) = ID;
                case(01)%y���򣬵ȼ�1
                    loadPf( x+1 ) = ID;
                    stifKf( x+1 ) = ID;
                case(11)%x��y����
                    loadPf( x ) = ID;
                    loadPf( x+1 ) = ID;
                    stifKf( x , x ) = ID;
                    stifKf( x+1 , x+1 ) = ID;
                otherwise
                    disp('DEBUG: LargeNumber!')
            end
        end
        stifKf(loadPf==ID,:) = [];%����
        stifKf(:,loadPf==ID) = [];%����
        loadPf(loadPf==ID) = [];%����
        loadP = loadPf;
        stifK = stifKf;
    otherwise
end
end



% xiezhuoyu
% mechanics_xzy@163.com
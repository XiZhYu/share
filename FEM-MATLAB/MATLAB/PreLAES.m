% 线性代数方程求解前处理
% 18/12/16
function [ loadP,stifK ] = PreLAES( loadPf,stifKf,fixN )
method = 1;
switch method
    %大数法
    case(1)
        LG = 1E+10;%大数
        for i = 1:1:size(fixN,1)
            x = (fixN(i,1)-1)*2 + 1;
            switch fixN(i,2)
                case(10)%x方向
                    loadPf( x ) = fixN(i,3)*LG;
                    stifKf( x,x ) = LG;
                case(01)%y方向，等价1
                    loadPf( x+1 ) = fixN(i,4)*LG;
                    stifKf( x+1,x+1 ) = LG;
                case(11)%x、y方向
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
    %划行法
    %破坏了矩阵原有数据、大小
    case(2)
        ID = 1E+30;%作为标记，方便后期整体删除
        for i = 1:1:size(fixN,1)
            x = (fixN(i,1)-1)*2 + 1;
            switch fixN(i,2)
                case(10)%x方向
                    loadPf( x ) = ID;
                    stifKf( x,x ) = ID;
                case(01)%y方向，等价1
                    loadPf( x+1 ) = ID;
                    stifKf( x+1 ) = ID;
                case(11)%x、y方向
                    loadPf( x ) = ID;
                    loadPf( x+1 ) = ID;
                    stifKf( x , x ) = ID;
                    stifKf( x+1 , x+1 ) = ID;
                otherwise
                    disp('DEBUG: LargeNumber!')
            end
        end
        stifKf(loadPf==ID,:) = [];%划行
        stifKf(:,loadPf==ID) = [];%划列
        loadPf(loadPf==ID) = [];%划行
        loadP = loadPf;
        stifK = stifKf;
    otherwise
end
end



% xiezhuoyu
% mechanics_xzy@163.com
%三角形3结点数据文件读入
%181215
function [ npoin,nelem,nvfix,ntype,nnode,nmats,nhamm,noutp,elemM,elemN, ...
    nodeC,fixN,mater,model,iplod,lodpt ] = InputFile( fileName )
fileID = fopen(fileName);
FEM_title = fgets(fileID);%读入第一行：说明
npoin = fscanf(fileID,'%d',[1,1]);%结点总数
nelem = fscanf(fileID,'%d',[1,1]);%单元总数
nvfix = fscanf(fileID,'%d',[1,1]);%约束数目
ntype = fscanf(fileID,'%d',[1,1]);%问题类型
nnode = fscanf(fileID,'%d',[1,1]);%单元结点数
nmats = fscanf(fileID,'%d',[1,1]);%材料数目
nhamm = fscanf(fileID,'%d',[1,1]);%Hammer积分阶次
noutp = fscanf(fileID,'%d',[1,1]);%结果输出控制参数
temp = fscanf(fileID,'%d',[nelem,2+nnode]);%读入单元信息
temp = reshape(temp(:),2+nnode,nelem);
temp = temp';
elemM = temp(:,2);%单元材料号
elemN = temp(:,3:1:2+nnode);%单元结点号
clear temp temp_nnode
temp = fscanf(fileID,'%d %f %f',[npoin,3]);%读入结点信息
temp = reshape(temp(:),3,npoin);
temp = temp';
nodeC = temp(:,2:1:3);%结点坐标
clear temp
temp = fscanf(fileID,'%d %d %f %f',[nvfix,4]);%读入约束信息
temp = reshape(temp(:),4,nvfix);
temp = temp';
fixN = temp;%约束信息
clear temp
temp = fscanf(fileID,'%d %f %f %f %f %d',[nmats,6]);%读入材料信息
temp = reshape(temp(:),6,nmats);
temp = temp';
mater = temp(:,1:1:5);%材料信息
model = temp(:,6);%单元模式号
clear temp
temp = fgets(fileID);%换行
clear temp
loads_title = fgets(fileID);%读入：载荷说明
iplod = fscanf(fileID,'%d %d %d',[1,3]);%读入载荷控制参数
temp = fscanf(fileID,'%d %f %f');%读入集中载荷
temp = reshape(temp(:),3,length(temp)/3);
temp = temp';
lodpt = temp;%集中载荷信息
clear temp
fclose(fileID);
end

% 数据格式
%   12 10  2 1 3 1 2 333
%      1     1     1     3     2
%      2     1     3     4     2
%      3     1     3     5     4
%      4     1     5     6     4
%      5     1     5     7     6
%      6     1     7     8     6
%      7     1     7     9     8
%      8     1     9    11     8
%      9     1     9    10    11
%     10     1    10    12    11
%      1     0.0000E+000     0.0000E+000
%      2     0.0000E+000     2.0000E+000
%      3     1.0000E+000     0.0000E+000
%      4     2.0000E+000     2.0000E+000
%      5     2.0000E+000     0.0000E+000
%      6     4.0000E+000     2.0000E+000
%      7     4.0000E+000     0.0000E+000
%      8     5.0000E+000     2.0000E+000
%      9     7.0000E+000     0.0000E+000
%     10     1.0000E+001     0.0000E+000
%     11     6.0000E+000     2.0000E+000
%     12     1.0000E+001     2.0000E+000
%      1        11     0.0000E+000     0.0000E+000
%      2        10     0.0000E+000     0.0000E+000
% 1
% 1 0.3 1.0 0.0 333
% 1 0 0
%      6     1.0000E+003     0.0000E+000
%     12    -1.0000E+003     0.0000E+000



% xiezhuoyu
% mechanics_xzy@163.com
%单元应力佳点，应力
%18/12/16
%单元结点应力
% 18/12/30
function [ strEP,strEN ] = Stres( asdis,DD,nelem,nhamm,elemM,elemN,nodeC,nnode,npoin )
posit = [   1/3 1/3  nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   ...
    ;   1/6 1/6  1/6 2/3 2/3 1/6 nan   nan   nan   nan   nan   nan   nan   nan   ...
    ;   1/3 1/3  1/5 1/5 1/5 3/5 3/5 1/5 nan   nan   nan   nan   nan   nan   ...
    ;   1/3 1/3  0.1012865073   0.1012865073    0.1012865073    0.7974269853    0.7974269853    0.1012865073 ...
    0.4701420641	0.4701420641	0.4701420641	0.0597158717	0.0597158717	0.4701420641];%方便走循环
%posit = [1/3 1/3 0 0 0 0 0 0 0 0 0 0 0 0;a a a b b a 0 0 0 0 0 0 0 0;1/3 1/3 a a a b b b 0 0 0 0 0 0;1/3 1/3 a a a b b a c c c d d c];
nPoiO = [1 3 4 0 7];%Hammer积分，各阶精度对应的积分点数目
strEP = zeros(3*nelem*nhamm,1);%3个应力分量、nelem个单元、nhamm个积分点
localP = [ 1 0 ; 0 1 ; 0 0 ];%结点局部坐标
strEN = zeros(3*npoin,1);%3个应力分量、nnode个结点
for i = 1:1:nelem%每个单元
    D = DD(:,:,elemM(i));%弹性矩阵
    for j = 1:1:nPoiO(nhamm)%每个积分点
        temp = (i-1)*nPoiO(nhamm)*3 + (nnode*(j-1)+1);%没什么意义，别太纠结。单元数*积分点数*应力分量数
        [ B,~ ] = StraB(i,posit(nhamm,(j-1)*2+1),posit(nhamm,(j-1)*2+2),nodeC,elemN,nnode);
        switch nnode
            case {3}
                strEP(temp:1:temp+2) = D*B ...
                    *[asdis((elemN(i,1)-1)*2+1) asdis((elemN(i,1)-1)*2+2) ...
                    asdis((elemN(i,2)-1)*2+1) asdis((elemN(i,2)-1)*2+2) ...
                    asdis((elemN(i,3)-1)*2+1) asdis((elemN(i,3)-1)*2+2) ]';
            case {6}
                strEP(temp:1:temp+2) = D*B ...
                    *[asdis((elemN(i,1)-1)*2+1) asdis((elemN(i,1)-1)*2+2) ...
                    asdis((elemN(i,2)-1)*2+1) asdis((elemN(i,2)-1)*2+2) ...
                    asdis((elemN(i,3)-1)*2+1) asdis((elemN(i,3)-1)*2+2) ...
                    asdis((elemN(i,4)-1)*2+1) asdis((elemN(i,4)-1)*2+2) ...
                    asdis((elemN(i,5)-1)*2+1) asdis((elemN(i,5)-1)*2+2) ...
                    asdis((elemN(i,6)-1)*2+1) asdis((elemN(i,6)-1)*2+2)]';
            otherwise
                disp('Stres1 ERROR!')
        end
        clear temp
    end
    for j = 1:1:nnode%每个结点
        [ B,~ ] = StraB(i,localP(j,1),localP(j,2),nodeC,elemN,nnode);
        switch nnode
            case {3}
                temp_streEN = D*B ...
                    *[asdis((elemN(i,1)-1)*2+1) asdis((elemN(i,1)-1)*2+2) ...
                    asdis((elemN(i,2)-1)*2+1) asdis((elemN(i,2)-1)*2+2) ...
                    asdis((elemN(i,3)-1)*2+1) asdis((elemN(i,3)-1)*2+2) ]';
                %绕点平均
                if(abs(strEN((elemN(i,j)-1)*nnode+1))<1E-6)
                    strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+3 ) =  temp_streEN;
                else
                    strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+3 ) = 0.5 * ( temp_streEN + strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+3 ) );
                end
            case {6}
                temp_streEN = D*B ...
                    *[asdis((elemN(i,1)-1)*2+1) asdis((elemN(i,1)-1)*2+2) ...
                    asdis((elemN(i,2)-1)*2+1) asdis((elemN(i,2)-1)*2+2) ...
                    asdis((elemN(i,3)-1)*2+1) asdis((elemN(i,3)-1)*2+2) ...
                    asdis((elemN(i,4)-1)*2+1) asdis((elemN(i,4)-1)*2+2) ...
                    asdis((elemN(i,5)-1)*2+1) asdis((elemN(i,5)-1)*2+2) ...
                    asdis((elemN(i,6)-1)*2+1) asdis((elemN(i,6)-1)*2+2)]';
                %绕点平均
                if(abs(strEN((elemN(i,j)-1)*nnode+1))<1E-6)
                    strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+6 ) = temp_streEN;
                else
                    strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+6 ) = 0.5 * ( temp_streEN + strEN( (elemN(i,j)-1)*nnode+1 :1: (elemN(i,j)-1)*nnode+6 ) );
                end
            otherwise
                disp('Stres2 ERROR!')
        end
    end
    clear temp
end
end


% xiezhuoyu
% mechanics_xzy@163.com
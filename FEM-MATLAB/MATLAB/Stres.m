%��ԪӦ���ѵ㣬Ӧ��
%18/12/16
%��Ԫ���Ӧ��
% 18/12/30
function [ strEP,strEN ] = Stres( asdis,DD,nelem,nhamm,elemM,elemN,nodeC,nnode,npoin )
posit = [   1/3 1/3  nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   ...
    ;   1/6 1/6  1/6 2/3 2/3 1/6 nan   nan   nan   nan   nan   nan   nan   nan   ...
    ;   1/3 1/3  1/5 1/5 1/5 3/5 3/5 1/5 nan   nan   nan   nan   nan   nan   ...
    ;   1/3 1/3  0.1012865073   0.1012865073    0.1012865073    0.7974269853    0.7974269853    0.1012865073 ...
    0.4701420641	0.4701420641	0.4701420641	0.0597158717	0.0597158717	0.4701420641];%������ѭ��
%posit = [1/3 1/3 0 0 0 0 0 0 0 0 0 0 0 0;a a a b b a 0 0 0 0 0 0 0 0;1/3 1/3 a a a b b b 0 0 0 0 0 0;1/3 1/3 a a a b b a c c c d d c];
nPoiO = [1 3 4 0 7];%Hammer���֣����׾��ȶ�Ӧ�Ļ��ֵ���Ŀ
strEP = zeros(3*nelem*nhamm,1);%3��Ӧ��������nelem����Ԫ��nhamm�����ֵ�
localP = [ 1 0 ; 0 1 ; 0 0 ];%���ֲ�����
strEN = zeros(3*npoin,1);%3��Ӧ��������nnode�����
for i = 1:1:nelem%ÿ����Ԫ
    D = DD(:,:,elemM(i));%���Ծ���
    for j = 1:1:nPoiO(nhamm)%ÿ�����ֵ�
        temp = (i-1)*nPoiO(nhamm)*3 + (nnode*(j-1)+1);%ûʲô���壬��̫���ᡣ��Ԫ��*���ֵ���*Ӧ��������
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
    for j = 1:1:nnode%ÿ�����
        [ B,~ ] = StraB(i,localP(j,1),localP(j,2),nodeC,elemN,nnode);
        switch nnode
            case {3}
                temp_streEN = D*B ...
                    *[asdis((elemN(i,1)-1)*2+1) asdis((elemN(i,1)-1)*2+2) ...
                    asdis((elemN(i,2)-1)*2+1) asdis((elemN(i,2)-1)*2+2) ...
                    asdis((elemN(i,3)-1)*2+1) asdis((elemN(i,3)-1)*2+2) ]';
                %�Ƶ�ƽ��
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
                %�Ƶ�ƽ��
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
function [ h,p ] = testNormality( x )
%TESTNORMALITY ��̬�ֲ�����
%   �˴���ʾ��ϸ˵��
    [h,p] = jbtest(x);
    disp(['p-value: ',num2str(p)])
    if h == 1
        disp('�ܾ���̬�ֲ���ԭ����')
    else
        disp('���ܾܾ���̬�ֲ���ԭ����')
    end

end


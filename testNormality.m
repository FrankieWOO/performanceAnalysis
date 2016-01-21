function [ h,p ] = testNormality( x )
%TESTNORMALITY 正态分布检验
%   此处显示详细说明
    [h,p] = jbtest(x);
    disp(['p-value: ',num2str(p)])
    if h == 1
        disp('拒绝正态分布的原假设')
    else
        disp('不能拒绝正态分布的原假设')
    end

end


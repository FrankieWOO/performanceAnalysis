function [ DD ] = downsideDeviation( x, target )
%UNTITLED downside deviation,i.e. downside risk
%   target: target return compared to x. Ŀ�������ʣ�ע�ⲻ���껯��MAR������x��Ҫ��ͬ����ʱ��߶��ϵ�
%
    %check if x is numeric 
    if ~isnumeric(x)
       error('x must be double. cannot accept other data types') 
    end
    % check whether target is numeric
    if ~isnumeric(target)
       error('target must be double. cannot accept other data types') 
    end
    % ���x�Ƿ�������
    if size(x,1) ~=1 && size(x,2) ~=1
        error('x must be a vector, cannot accept matrix in this version of code')
    end
    % ǰ���Ѿ�������Ƿ�������
    n = length(x);
    d = max(target - x,0);
    DD = sqrt(sum( d.^2 ))/n;
    
end


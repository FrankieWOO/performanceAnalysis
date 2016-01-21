function [ DD ] = downsideDeviation( x, target )
%UNTITLED downside deviation,i.e. downside risk
%   target: target return compared to x. 目标收益率，注意不是年化的MAR，它和x需要是同样的时间尺度上的
%
    %check if x is numeric 
    if ~isnumeric(x)
       error('x must be double. cannot accept other data types') 
    end
    % check whether target is numeric
    if ~isnumeric(target)
       error('target must be double. cannot accept other data types') 
    end
    % 检查x是否是向量
    if size(x,1) ~=1 && size(x,2) ~=1
        error('x must be a vector, cannot accept matrix in this version of code')
    end
    % 前面已经检查了是否是向量
    n = length(x);
    d = max(target - x,0);
    DD = sqrt(sum( d.^2 ))/n;
    
end


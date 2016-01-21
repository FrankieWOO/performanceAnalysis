function [ B ] = cleanNanRow( A , var )
%cleanNanRow 清理table中var为NaN的行,var为variableName
%   A.(var)的用法是table的动态变量名引用
    if iscellstr(A.(var) ) || iscell(A.(var))
        x = strcmp( A.(var), '' ) ;
    else
        x = isnan( A.(var) ) ;
    end
    B=A(~x,:);
    disp(['clean ',num2str(sum(x)),' Na rows'])
end


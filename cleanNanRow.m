function [ B ] = cleanNanRow( A , var )
%cleanNanRow ����table��varΪNaN����,varΪvariableName
%   A.(var)���÷���table�Ķ�̬����������
    if iscellstr(A.(var) ) || iscell(A.(var))
        x = strcmp( A.(var), '' ) ;
    else
        x = isnan( A.(var) ) ;
    end
    B=A(~x,:);
    disp(['clean ',num2str(sum(x)),' Na rows'])
end


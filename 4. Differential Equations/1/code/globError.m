function [error] = globError(Y,ynList)
[m,timeLen] = size(ynList);
error = zeros(1,timeLen);
for i = 1:timeLen
    error(i) = norm(Y(:,i)-ynList(:,i));
end
end


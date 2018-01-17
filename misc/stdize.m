

function [X,stobj] = stdize(X)

stobj.mu = mean(X);
stobj.st = std(X);

idx = find(stobj.st>0);
% keyboard
X(:,idx) = bsxfun(@minus,X(:,idx),stobj.mu(idx));
X(:,idx) = bsxfun(@rdivide,X(:,idx),stobj.st(idx));

% for i=1:size(X,2)
%     if stobj.st(i)~=0
%         X(:,i) = (X(:,i)-stobj.mu(i))/stobj.st(i);
%     end
% end
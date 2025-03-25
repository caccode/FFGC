function [C, F, y] = BalancedKM_plus(X, ratio)
class_num = 2;
n = size(X,2);
if nargin < 2
    ratio = 0.5;
end
%   StartInd = randsrc(n,1,1:class_num);
C=X(:,1+round(rand*(n-1)));
StartInd=ones(1,n);
for i=2:class_num
    D=X-C(:,StartInd);
    D=cumsum(sqrt(dot(D,D,1)));
    if D(end)==0
        C(:,i:class_num)=X(:,ones(1,class_num-i+1));
    else
    C(:,i)=X(:,find(rand<(D/D(end)),1));
    end
    [~,StartInd]=max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    StartInd=StartInd';
end
InitF = TransformL(StartInd, class_num);

if ratio > 0.5
    error('ratio should not larger than 0.5');
end
if ratio < 0
    ratio = 0;
end

a = floor(n*ratio);

F = InitF;
last = InitF(:,1);
F_zeros = zeros(n,class_num);
for iter = 1:100
    C = X*F*inv(F'*F+eps*eye(2));
    F = F_zeros;
    Q = L2_distance_1(X,C);
    q = Q(:,1)-Q(:,2);
    [temp, idx] = sort(q);
    cp=a;
    F(idx(1:cp),1) = 1;
    F(:,2) = 1-F(:,1);
    if F(:,1)==last(:)
        break;
    end
    last = F(:,1);
end

[~, y] = max(F,[],2);
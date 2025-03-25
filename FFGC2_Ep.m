function [F,objs,runtime] = FFGC2_Ep(F,B,c,gamma)
% FFGC2 Fast Fuzzy Graph cluertering _ reweighted optimization
% This code is implemented by Qianyao Qiang
% B in n by m
% F in n by c
% s in c by 1
NITR=100;
NITR_RW=100;
[n,~]=size(F);
s=zeros(c,1);BF2=zeros(c,1);FW=zeros(c,1);

% pre calculate
sumF=sum(F);
BF=B'*F; % m by c
for j=1:c
    BF2(j)=(BF(:,j)'*BF(:,j))^(1/2);
end

tic
for iter=1:NITR
    % update s
    for j=1:c
        s(j)=BF2(j)/sumF(j);
    end
    
    % update F    % reweighted
    s2=s.^2;% c by 1
    for iter2=1:NITR_RW
        % update weight
        BF2_=diag(1./BF2);% c by c
        W=B*(BF*BF2_);% W in n by c
        
        % update F
        E=2*W*diag(s);% n by c
        for i=1:n
            hnew=(E(i,:)-s2')/(2*gamma);
            F(i,:)=EProjSimplex_new(hnew);
        end
        
        BF=B'*F;
        sumF=sum(F);
        for j=1:c
            BF2(j)=(BF(:,j)'*BF(:,j))^(1/2);
            FW(j)=F(:,j)'*W(:,j);
        end
        
        objw(iter2)=2*s'*FW-(s.^2)'*sumF'-gamma*norm(F,'fro')^2;
        if iter2>2 && abs((objw(iter2)-objw(iter2-1))/objw(iter2))<1e-10
            break;
        end
        if iter2>30 && sum(abs(objw(iter2-9:iter2-5)-objw(iter2-5+1:iter2)))<1e-10
            break;
        end
    end
    
    objs(iter)= 2*s'*BF2-(s.^2)'*sumF'-gamma*norm(F,'fro')^2;
    if iter>2 && abs((objs(iter)-objs(iter-1))/objs(iter))<1e-10
        break;
    end
    if iter>30 && sum(abs(objs(iter-9:iter-5)-objs(iter-5+1:iter)))<1e-10
        break;
    end
end
runtime=toc;
end


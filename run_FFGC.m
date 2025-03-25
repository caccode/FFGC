function [result,runtime,obj] = run_FFGC(ds,numAnchor,gamma,s)
% Run Fast Fuzzy Graph cluertering
load (ds);
X = full(X);  % num*dim
if min(Y)==1
    gnd = Y;
else
    gnd = Y + 1-min(Y);
end
c = max(gnd);
numNearestAnchor = 5;
X =(double(X));  %n*d
n = size(X,1);
%% normalization
    % centered
%      for  j = 1:n
%         X(j,:) = ( X(j,:) - mean( X(j,:) ) ) ;
%     end
    % max norm
%      a = max(X(:));
%      X = double(X./a);
    % standard norm
%     X = (X-repmat(mean(X),[n,1]))/repmat(std(X),[n,1]);
    % max and min norm
%     X=( X-repmat(min(X),n,1) )./repmat(max(X)-min(X),n,1);
%     X(isnan(X))=1;
%% construct B
[~,locAnchor] = hKM(X',[1:n],numAnchor,1); % here we use BKHK algorithm with k-means repleacing with k-means++.
Z = ConstructA_NP(X',(locAnchor),numNearestAnchor); 
sumZ = sum(Z);
sqrtZ = sumZ.^(-0.5);
B = (Z)*(diag(sqrtZ));
%% init F
% r=1.5; [prediction,v,U,obj_Fcm] = fcm(c, X, r); F=U';
F=rand(n,c); F=F./repmat(sum(F,2),[1,c]);
if s==1
    [F,obj,runtime] = FFGC1(F,B,c,gamma);
elseif s==2
    [F,obj,runtime] = FFGC2_Ep(F,B,c,gamma);
end
[~,F]=max(F,[],2);
result = ClusteringMeasure(gnd,F);
end


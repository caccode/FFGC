clear all
close all
clc
%--------------------------------------------------------------------------------
addpath('datasets','functions','measure');
ds = 'isolet_uni.mat';
s = 2; % 1;
fprintf('Result on %s \r\n', ds);
for numAnchor=[6,7,8,9,10,11]
    %% Fast Fuzzy Graph Clustering
    for gamma=[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2]
        ACC_ours = []; NMI_ours = []; Purity_ours = []; T_ours = [];
        for cnt=1:20
            [result,runtime,obj] = run_FFGC(ds,numAnchor,gamma,s);
            ACC_ours(cnt,:)=result(1); NMI_ours(cnt,:)=result(2); Purity_ours(cnt,:)=result(3);T_ours(cnt,:)=runtime;
            fprintf('ours ACC = %0.4f NMI = %0.4f Purity = %0.4f RunTime= %f numAnchor=%f gamma=%f \r\n',result,runtime,numAnchor,gamma);
        end
        fprintf('ours1_meanresult ACC = %0.4f NMI = %0.4f  Purity = %0.4f RunTime= %0.4f  numAnchor=%f gamma=%f \r\n',mean(ACC_ours),mean(NMI_ours),mean(Purity_ours),mean(T_ours),numAnchor,gamma);
        fprintf('ours1_stdresult ACC = %0.4f NMI = %0.4f  Purity = %0.4f RunTime= %0.4f   numAnchor=%f gamma=%f \r\n',std(ACC_ours),std(NMI_ours),std(Purity_ours),std(T_ours),numAnchor,gamma);
    end

end


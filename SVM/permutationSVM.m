function [classAcc, accs, p95threshold] = permutationSVM(svmStruct, betas, permLabels, runs)
% betas: trials x features matrix
% permLabels: permuted labels following permuteLabels.m for trials
% runs: run labels
% trials: presentation number. (permute the 1s, permute the 2s, etc etc)
% PLEASE PLEASE PLEASE HAVE A BALANCED NUMBER OF TRIALS ACROSS RUNS
% runs a proper permutation test for cross validation


% for testing
% betas=[1*rand(50, 300); 2*rand(50,300); 1*rand(50, 300); 2*rand(50,300)];
% labels=[ones(50,1); 2*ones(50,1);ones(50,1); 2*ones(50,1)];
% runs=[ones(100,1); 2*ones(100,1)];
% % trials=[[1:50]'; [1:50]'; 50+[1:50]'; 50+[1:50]'];

nPerm=size(permLabels,2);
nRuns=length(unique(runs));
opts=['-s 1 -t 0'];

try  
    %% parallelize across subjects
    
    accs=zeros(nRuns, nPerm);
    permI=1;
    %% run the classifications
    for permI=1:nPerm
        labelIdx=permLabels(:,permI);
        for runI = 1:nRuns
            % for cross validation
%             trainPatterns=betas(runs~=runI,:);
%             trainLabels=labelIdx(runs~=runI,:);
%             svmStruct=libsvmtrain(trainLabels, trainPatterns, opts);
            testPatterns=betas(runs==runI,:);
            testLabels=labelIdx(runs==runI,:);
            [~, acc, ~] = libsvmpredict(testLabels, testPatterns, svmStruct);
            accs(runI, permI)=acc(1);
        end
    end
    
    %% get stats
    % this should give the distribution/percentile
    p95threshold=quantile(mean(accs), .95);
    classAcc=mean(accs(:,1));
%     fprintf('THE CLASSIFICATION ACCURACY IS DUR... %d \n', classAcc)
%     fprintf('THE .05 CUTOFF IS DUR... %d \n', p95threshold)
catch err
    save('permutationError.mat')
    rethrow(err)
end
end

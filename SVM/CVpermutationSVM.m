function [result] = CVpermutationSVM(betas, permLabels, runs, varargin)
% function [] = CVpermutationSVM(betas, permLabels, runs, varargin)
% betas: trials x features matrix
% permLabels: permuted labels following permuteLabels.m for trials
% runs: run labels
% runs a proper permutation test for cross validation (leave one run out)


%% parse input
p=inputParser;
defaultOpts=['-q -s 1 t 0'];
% validFinishes = {'-s 1 t 0','-s 1 -t 2', '-s 0 -t 0', '-s 0 -t 2'};
% checkOpts = @(x) any(validatestring(x,validFinishes));
% defaultKernel=0;
addRequired(p,'betas', @isnumeric);
addRequired(p,'permLabels', @isnumeric);
addRequired(p,'runs', @isnumeric);
addParameter(p,'opts', defaultOpts, @ischar);
% addParameter(p,'svmKernel', defaultKernel, @isnumeric);

p.KeepUnmatched=true;
parse(p, betas, permLabels, runs, varargin{:});
% applyZ=@(betas, mu, sigma) (bsxfun(@times,bsxfun(@minus,betas, mu),1./sigma));

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end
if ~isempty(p.UsingDefaults)
    disp('Using defaults: ')
    disp(p.UsingDefaults)
end

%% initialize parameters
opts=p.Results.opts;
% svmKernel=p.Results.svmKernel;
nPerm=size(permLabels,2);
includeIdx=[permLabels(:,1)>0];
betas=betas(includeIdx,:); % exclude trials with 0 as label
runs=runs(includeIdx,:);
nRuns=length(unique(runs));

try
    % initialize and zscore
    result.acc=zeros(nRuns, nPerm);
    result.mse=zeros(nRuns, nPerm);
    result.scc=zeros(nRuns, nPerm);
    origLabels=permLabels(:,1);
    
    for r = 1:nRuns
        betas(runs==r,:)=zscore(betas(runs==r,:));
    end
    
    %% run the classifications
    if nRuns==1
        % if a all in a single run, just do 2-fold cross validation
        for permI=1:nPerm
            labelIdx=permLabels(:,permI);
            trainPatterns=zscore(betas);
            trainLabels=labelIdx;
            result.acc(1,permI)=libsvmtrain(trainLabels, trainPatterns, [opts ' -v 2']);
            disp('Running 2 fold using libsvm default. Are you sure?')
        end
    else
        for runI = 1:nRuns
            %% set up training and testing set for this iteration
            % soft normalization by chunks
            trainLabels=origLabels(runs~=runI,:);
            trainPatterns=betas(runs~=runI,:);
            testPatterns=betas(runs==runI,:);
            % not doing optimize SVM
            svmStruct=libsvmtrain(trainLabels, trainPatterns, opts);
            % permutation: test on random labels
            for permI=1:nPerm
                labelIdx=permLabels(:,permI);
                testLabels=labelIdx(runs==runI,:); % classify on original labels, but shuffle dinput
                [~, acc, ~] = libsvmpredict(testLabels, testPatterns, svmStruct);
                result.acc(runI, permI)=acc(1);
                result.mse(runI,permI)=acc(2);
                result.scc(runI,permI)=acc(3);
            end
            
        end
    end
    
    %% get stats
    % this should give the distribution/percentile
    result.p95threshold=quantile(mean(result.acc,1), .95);
    result.classAcc=mean(result.acc(:,1));
catch err
    save('permutationError.mat')
    rethrow(err)
end
end


function [classAcc, accs, p95threshold] = CVpermutationSVM(betas, permLabels, runs, varargin)
% function [] = CVpermutationSVM(betas, permLabels, runs, varargin)
% betas: trials x features matrix
% permLabels: permuted labels following permuteLabels.m for trials
% runs: run labels
% PLEASE PLEASE PLEASE HAVE A BALANCED NUMBER OF TRIALS ACROSS RUNS
% runs a proper permutation test for cross validation (leave one run out)


%% parse input
p=inputParser;
defaultOpts=['-s 1 t 0'];
% validFinishes = {'-s 1 t 0','-s 1 -t 2', '-s 0 -t 0', '-s 0 -t 2'};
% checkOpts = @(x) any(validatestring(x,validFinishes));
defaultKernel=0;
addRequired(p,'betas', @isnumeric);
addRequired(p,'permLabels', @isnumeric);
addRequired(p,'runs', @isnumeric);
addParameter(p,'opts', defaultOpts, @ischar);
addParameter(p,'svmKernel', defaultKernel, @isnumeric);

p.KeepUnmatched=true;
parse(p, betas, permLabels, runs, varargin{:});
applyZ=@(betas, mu, sigma) (bsxfun(@times,bsxfun(@minus,betas, mu),1./sigma));

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
svmKernel=p.Results.svmKernel;
nPerm=size(permLabels,2);
includeIdx=[permLabels(:,1)>0];
betas=betas(includeIdx,:); % exclude trials with 0 as label
runs=runs(includeIdx,:);
nRuns=length(unique(runs));

try
    accs=zeros(nRuns, nPerm);
    permI=1;
    origLabels=permLabels(:,1);
    %% run the classifications    
    if nRuns==1
        % if a all in a single run, just do 2-fold cross validation
        for permI=1:nPerm
            labelIdx=permLabels(:,permI);
            trainPatterns=zscore(betas);
            trainLabels=labelIdx;
            accs(1,permI)=libsvmtrain(trainLabels, trainPatterns, [opts ' -v 2']);
            disp('Running 2 fold using libsvm default. Are you sure?')
        end
    else
        for runI = 1:nRuns
%% set up training and testing set for this iteration
% zscore each training set, and use that mean and SD for the

% testing set.
            [trainPatterns, timeMu, timeSigma]=zscore(betas(runs~=runI,:));
            trainLabels=origLabels(runs~=runI,:);
            testPatterns=applyZ(betas(runs==runI,:), timeMu, timeSigma);

            % optimize SVM
            opt_params=optimizeSVM(trainLabels, trainPatterns, svmKernel);
            if svmKernel==2
                opts=['-q -s 0 -t ' num2str(svmKernel) ' -c ' num2str(opt_params.best_C) ' -g ' num2str(opt_params.best_gamma)];
            else
                opts=['-q -s 1 -t ' num2str(svmKernel) ' -n ' num2str(opt_params.best_C) ];
            end
            svmStruct=libsvmtrain(trainLabels, trainPatterns, opts);
            
            % permutation: test on random labels
            for permI=1:nPerm
                labelIdx=permLabels(:,permI);
                testLabels=labelIdx(runs==runI,:); % classify on original labels, but shuffle dinput
                [~, acc, ~] = libsvmpredict(testLabels, testPatterns, svmStruct);
                accs(runI, permI)=acc(1);
            end
            
        end
    end
    
    %% get stats
    % this should give the distribution/percentile
    p95threshold=quantile(mean(accs,1), .95);
    classAcc=mean(accs(:,1));
catch err
    save('permutationError.mat')
    rethrow(err)
end
end
% this part is wrong-- normalizing it by space then time removes important
% differences.
%             [trainPatterns, spaceMu, spaceSigma]=zscore(betas(runs~=runI,:),0,2);
%             [trainPatterns, timeMu, timeSigma]=zscore(trainPatterns);
%             zscore across time. NOT space.
%             trainPatterns=betas(runs~=runI,:);
%             testPatterns=zscore(betas(runs==runI,:), 0, 2);
%             testPatterns=applyZ(testPatterns, timeMu, timeSigma);
%             testPatterns=betas(runs==runI,:);

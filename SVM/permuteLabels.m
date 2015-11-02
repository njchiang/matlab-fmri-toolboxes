function [permLabels] = permuteLabels(labels, runs, nPerm)
% labels: labels for trials
% runs: run labels
% trials: presentation number. (permute the 1s, permute the 2s, etc etc)
% PLEASE PLEASE PLEASE HAVE A BALANCED NUMBER OF TRIALS ACROSS RUNS
% runs a proper permutation test for cross validation
% permuted_labels: one label in each column
%  later: add capability to use presentation order


%% permute the labels (using runs and trials)
% ONE PERMUTATION, CONSISTENT ACROSS SUBJECTS
% for now, permutes all trials within a run. will eventually code super
% specific ones that keep presentation constant
permLabels=zeros(length(labels), nPerm);
% first one is the original
origOrder=[1:length(labels)]';
runOrig=[runs origOrder];
permLabels(:,1)=labels;
for i = 2:nPerm
    %     permute label for given trial number in given run
    labelIdx=randperm(length(labels))';
    permutedLabels=sortrows(runOrig(labelIdx,:),1);
    theseLabels=labels(permutedLabels(:,2));
    permLabels(:,i)=theseLabels;
end



% defineUserOptions;

toolboxRoot = '/space/raid5/data/monti/Analysis/LanguageMVPA/RSA/code'; addpath(genpath(toolboxRoot));
% userOptions = defineUserOptions_Syntax();
userOptions = defineUserOptions_LSemantics();
gotoDir(userOptions.rootPath);
analysisType='SVM';
userOptions.analysisType=analysisType;
%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%
% load data and masks
% fullBrainVols = fMRIDataPreparation('SPM', userOptions);

%betaCorrespondence_SVM in the order of runs.
fullBrainVols=fMRIDataPreparation(betaCorrespondence_(), userOptions);
binaryMasks_nS = fMRIMaskPreparation(userOptions);
responsePatterns = fMRIDataMasking(fullBrainVols, binaryMasks_nS, betaCorrespondence_SVM(), userOptions);
clear fullBrainVols binaryMasks_nS

%define runs, labels

%% permutation testing
%generate permutations
models=makeLabels_Semantics();
runs= %something related to betaCorrespondence
for i = 1:length(models)
    
    modelMat{i}=permuteLabels(models(i).label, runs, nPerm);
end


parpool open
parfor j=1:length(userOptions.subjectNames)
currSub=userOptions.subjectNames{j}
for m =1:length(userOptions.maskNames)
    currMask=userOptions.maskNames{m}
    thesePatterns=responsePatterns.(currMask).(currSub)';
    if size(thesePatterns,1)~=size(modelMat,1)
        thesePatterns=thesePatterns'
    end
    singleSubjectPattern=zscore(thesePatterns, 0, 2); % make sure this is right
    permuteSVM(singleSubjectPattern, modelMat, runs);
end
end


searchlightOptions.monitor = false;
searchlightOptions.fisher = true;

nCond = searchlightOptions.nConditions;


%PERMUTE LABELS
[smm_rs, smm_ps, n, searchlightRDMs] = searchlightMapping_SVM(fullBrainVolumes, modelMAT, mask, runs, userOptions, searchlightOptions);
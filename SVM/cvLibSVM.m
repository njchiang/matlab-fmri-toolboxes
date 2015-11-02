function [results]=cvLibSVM(L,X,varargin)
% function to perform SVM on data X over labels L, doing an exhaustive N-1
% cross validation.
%
% syntax [svm_w]=cvLibSVM(X,L)
%       INPUT:
%       X is a data matrix with as many rows as observations and columns as
%       factors. NOTE: it might be worth using the command as
%       mysvam(zscore(X),L).
%       L is a label vector with as many rows as X.
%
%       OPTIONAL:
%           'nstochastic' is the number of iterations you want to run (if not
%           exhaustive
%           'opts' is a string with parameters for libSVM. default is -s 1
%           -t 0 (linear classification) NOTE: DO NOT specify cost. other
%           options below.
%           'display_coords' is a coordinate matrix for plotting the result importance map in
%           3D using scatter 3D
%           'margin' is the cost parameter, default is 1
%
%       OUTPUT:
%       results is a struct containing prediction, accuracy, and classifier
%       output in that order.
%
%        SVM_Options:
% -s svm_type : set type of SVM (default 0)
% 	0 -- C-SVC		(multi-class classification)
% 	1 -- nu-SVC		(multi-class classification)
% 	2 -- one-class SVM
% 	3 -- epsilon-SVR	(regression)
% 	4 -- nu-SVR		(regression)
% -t kernel_type : set type of kernel function (default 2)
% 	0 -- linear: u'*v
% 	1 -- polynomial: (gamma*u'*v + coef0)^degree
% 	2 -- radial basis function: exp(-gamma*|u-v|^2)
% 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
% 	4 -- precomputed kernel (kernel values in training_set_file)
% -d degree : set degree in kernel function (default 3)
% -g gamma : set gamma in kernel function (default 1/num_features)
% -r coef0 : set coef0 in kernel function (default 0)
% -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
% -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
% -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
% -m cachesize : set cache memory size in MB (default 100)
% -e epsilon : set tolerance of termination criterion (default 0.001)
% -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
% -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
% -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
% -v n: n-fold cross validation mode
% -q : quiet mode (no outputs)


%% parse inputs
defaults = struct('nstochastic', 0, 'opts', ['-s 1 -t 0'],'display_coords',[], 'margin', 1, 'nu', .2);
if nargin >=3;
    
    options = libparse_input( defaults, varargin );
    
else
    options = defaults;
end
opts = options.opts;
nStochastic = options.nstochastic;
V = options.display_coords;
C = options.margin;
nu = options.nu;
%% initialize
opts = [ opts ' -q -c ' num2str(C) ' -n ' num2str(nu)];
tp=0;tn=0;fp=0;fn=0;w=0; % set various counters to zero
if nStochastic==0
    combs = generate_iterations(L);
else
    try
        combs = generate_iterations(L,nStochastic);
    catch err
        disp(['More iterations requested than combinations, running max'])
        combs=generate_iterations(L);
    end
end
results.options = struct('opts', opts, 'Display_coords', V, 'nStochastic', nStochastic);

pred_c = [];
acc_c = [];
est_c = [];


disp(['**********************************'])
disp([num2str(size(combs,1)) ' iterations to be performed'])
disp(['Using opts: ' opts]);
disp(['**********************************'])
numIterations = size(combs,1);


for i=1:numIterations
    %% initalize within iteration
    extract_i = combs(i,:);
    train_X=X;       % Put the full set of data (to then hold 1 row out)
    train_L=L;       % Put all lables here (to then hold 1 row out)
    train_X(extract_i,:)=[]; % Pull out the ith row
    train_L(extract_i)=[];   % Pull out the ith row
    test_X=X(extract_i,:);   % Put the ith row in the testing set
    test_L=L(extract_i); % Put the ith label in the testing set
    %% train and test
    svmStruct=libsvmtrain(train_L,train_X,opts);
    
    [pred acc est] = libsvmpredict(test_L, test_X, svmStruct);
    est_c = [est_c est];
    acc_c = [acc_c acc];
    pred_c = [pred_c pred];
    %     if (pred(i) == 0) && (test_L == 0)
    %         tn=tn+1; %correctly said "negative"
    %         w=w+1;   %Up the counter to put the correct classification weights
    %         %        [svm_w(w,:),svm_b(w)]=svmcoefs(svmStruct);
    %     elseif (pred(i) == 1) && (test_L == 1)
    %         tp=tp+1; %correctly said "positive"
    %         w=w+1;   %Up the counter to put the correct classification weights
    %         %        [svm_w(w,:),svm_b(w)]=svmcoefs(svmStruct);
    %     elseif (pred(i) == 0) && (test_L == 1)
    %         fn=fn+1; %incorrectly said "negative"
    %     elseif (pred(i) == 1) && (test_L == 0)
    %         fp=fp+1; %incorrectly said "positive"
    %     end
end

%% store results
results.pred=pred_c;
results.acc = acc_c;
results.est = est_c;
results.svmStruct = libsvmtrain(L,X, opts);
results.accuracy = mean(results.acc(1,:));

disp(['Total Accuracy: ' num2str(results.accuracy)]);
%% display
if size(options.display_coords,2)== size(X,2)
    visualize_binary_svm(X,results.options.Display_coords, results.svmStruct, L);
end


%% legacy stats
%
% 3. Print out the confusion matrix
%     tot=size(L,1);
%     p_tot=sum(L==1);
%     n_tot=tot-p_tot;
%     % see http://en.wikipedia.org/wiki/Sensitivity_and_specificity
%     disp(' ')
%     disp('*******************************')
%     disp(['Classification performance: ']);% char(roi(r))])
%     disp('*******************************')
%     disp(['Classification Accuracy: ' num2str((tn+tp)/size(L,1))]);
%     disp(['Number of correct classifications: ' num2str(w/size(X,1))]);
%     disp(['  Sensitivity (true ''+'' rate/recall): ' num2str((tp)/(tp+fn))]); %True positive rate
%     disp(['  Specificity (true ''-'' rate): ' num2str((tn)/(tn+fp))]); %True Negative Rate
%     disp(['  Precision   (''+'' predic val):' num2str((tp)/(tp+fp))]); %Precision or positive predictive value
%     disp(' ' )
%     disp(['                        True:'])
%     disp(['                   -             +'])
%     disp(['------------------------------------------'])
%     disp([' Prediction - | ' num2str(tn/n_tot) '    ' num2str(fn/p_tot)])
%     disp(['            + | ' num2str(fp/n_tot) '    ' num2str(tp/p_tot)])
%     disp(['------------------------------------------'])
%         results.confusion_matrix = [tn/n_tot fn/p_tot; fp/n_tot tp/p_tot];
%             results.precision = (tp)/(tp+fp);
%     results.specificity = (tn)/(tn+fp);
%     results.sensitivity = (tp)/(tp+fn);
%     results.accuracy = (tn+tp)/size(L,1);



end

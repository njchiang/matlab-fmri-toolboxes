function [optimized_params] = optimizeSVM(trainLabels, trainPatterns, opts)
% UNDER CONSTRUCTION...
%grid of parameters
try
    svmKernel = opts(find(opts=='t')+2);
catch
    disp('could not find svmKernel, using 0')
end
folds = 4;
if svmKernel==2
[C,gamma] = meshgrid(-5:2:15, -15:2:3); 
%# grid search, and cross-validation
cv_acc = zeros(numel(C),1);
d= svmKernel;
for i=1:numel(C)
    cv_acc(i) = libsvmtrain(trainLabels,trainPatterns, ...
        sprintf('-c %f -g %f -v %d -t %d', 2^C(i), 2^gamma(i), folds,d));
end
%# pair (C,gamma) with best accuracy
[~,idx] = max(cv_acc);
% %# contour plot of paramter selection
% contour(C, gamma, reshape(cv_acc,size(C))), colorbar
% hold on;
% text(C(idx), gamma(idx), sprintf('Acc = %.2f %%',cv_acc(idx)), ...
%     'HorizontalAlign','left', 'VerticalAlign','top')
% hold off
% xlabel('log_2(C)'), ylabel('log_2(\gamma)'), title('Cross-Validation Accuracy')
optimized_params.best_C = 2^C(idx);
optimized_params.best_gamma = 2^gamma(idx); %# ...
optimized_params.best_acc=max(cv_acc);
else
[C] = [0.01:.05:.99]; 
%# grid search, and cross-validation
cv_acc = zeros(numel(C),1);
d= svmKernel;
for i=1:numel(C)
    cv_acc(i) = libsvmtrain(trainLabels,trainPatterns, ...
        sprintf('-n %f -v %d -t %d', C(i), folds,d));
end
%# pair (C,gamma) with best accuracy
[~,idx] = max(cv_acc);
%# contour plot of paramter selection
% plot(C, cv_acc)
% hold on;
% hold off
% xlabel('Nu'), ylabel('P'), title('Cross-Validation Accuracy')
optimized_params.best_C = C(idx);
optimized_params.best_acc=max(cv_acc);

optimized_params = '';
end





end

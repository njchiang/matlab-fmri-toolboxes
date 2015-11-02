function [ labels_perm ] = generate_iterations( labels, n )
%Takes label vector, orders it from 1:numClasses, and permutes it for
%exhaustive or stochastic cross validation of n samples
[~,~,L] = unique(labels);
perms = [];
for class = 1:max(L)
    perms(:,class) = find(L==class);
end
command_c = [];

for c = 1:max(L)-1
    command_c=[command_c 'perms(:,' mat2str(c) '), '];
end
command_c = [command_c 'perms(:,' mat2str(max(L)) ')'];
comb_stimuli = eval(['allcomb(' command_c ');']);
% exhaustive
numIterations = size(comb_stimuli,1);
labels_perm = comb_stimuli;


%stochastic search, exhaustive too long
if nargin==2
    randIterations = randsample(1:size(comb_stimuli,1), n);
    labels_perm = comb_stimuli(randIterations,:);
end

end


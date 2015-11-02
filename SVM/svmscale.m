function [scaled_data, scaling_factors] = svmscale(train_patterns, varargin)
% train patterns must be in correct format, numConditions x numFeatures
% each feature will be scaled to [-1 1]
% apply scaling_factors by multiplying train_patterns * scaling_factors;
if nargin <=4
    if nargin==4
        disp(['Applying previously generated scaling factors']);
        lbound=varargin{1};
        ubound=varargin{2};
        scaling_factors=varargin{3};
        for j = 1:size(train_patterns,2)
            x=train_patterns(:,j);
            scaled_data(:,j)=(x-min(x))*scaling_factors(j,j) + lbound;
        end
    else
        if nargin==1
            ubound=1;
            lbound=0;
            disp(['Lower bound and upper bound not specified, defaulting to [0 1]']);
            
        elseif nargin==3
            % x'=(upper-lower)*(x-min(x)/(max(x)-min(x)) + a;
            % matrix solution
            lbound=varargin{1};
            ubound=varargin{2};
            disp(['setting lower bound to ' num2str(lbound) ' and upper bound to ' num2str(ubound)]);
        end
        
        scaling_factors=(ubound-lbound)./(max(train_patterns)-min(train_patterns));
        scaling_factors=diag(scaling_factors);
        % scaled_data=(train_patterns-repmat(min(train_patterns), [size(train_patterns,1) 1]) * scaling_factors) + lower;
        for j = 1:size(train_patterns,2)
            x=train_patterns(:,j);
            %     scale_factor=(ubound-lbound)/(max(x)-min(x));
            scaled_data(:,j)=(x-min(x))*scaling_factors(j,j) + lbound;
        end
    end
    
else
    disp(['WHAT ARE YOU DOING?!?!'])
end

end
function [ options ] = libparse_input( defaults, new_args )
% function to parse inputs for SVM
% 
%# define defaults at the beginning of the code so that you do not need to
%# scroll way down in case you want to change something or if the help is
%# incomplete
% 
options = defaults;
% will need to add correction method for c
%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(new_args);
if round(nArgs/2)~=nArgs/2
   error('cvLibSVM needs propertyName/propertyValue pairs')
end

for pair = reshape(new_args,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive

   if any(strmatch(inpName,optionNames))
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end
end


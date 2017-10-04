function [coords_struct] = fmri_combine_coords(sessionFile,run_idx)

%
% USAGE: [coords_struct] = fmri_combine_coords(sessionFile,run_idx)
%
%   Combine the coords of several runs to a single coords, which is the
%   index of brain voxels that exist for all runs. The coords are assumed
%   stored in ${pls_data_path}/${prefix_datamat}_run?.mat files.
%
%   NOTE: all coords must come from the realigned datamats for the same 
%         subject. 
%   
%   Input: 
%      sessionFile:  the name of file that contains the session information.
%      run_idx:    the indices of runs that will be used to determine the 
%      (optional)  common coords.  If run_idx is not specified all coords 
%                  are used.
%
%   Output:
%      coords_struct - structure with the following 2 fields
%          common: the common coords 
%          coords_idx: a cell array contains the indices of the individual 
%                      coords for a run that can be used to generate the 
%	 	       common coords.
%
%   Output File: (** only generate if no output argument is provided **)
%      ${pls_data_path}/${prefix_datamat}_common_coords.mat:  
%           coords_info - (same as the above coords_struct)
%           run_idx: the indices of the runs from which the data are used 
%		     to compute the common coords
%
%      NOTE:  the new datamat for the common coords can be generated by
%	      new_datamat = datamat(coords_idx_run?);
%
%   EXAMPLE:
%       coords_struct = fmri_combine_coords('PLSsession',[1:3]);
%

   progress_hdl = ShowProgress('initialize');

   load(sessionFile);		% load the session information

   pls_data_path = session_info.pls_data_path;

   num_runs = session_info.num_runs;
   if ~exist('run_idx','var') | isempty(run_idx),
      run_idx = [1:num_runs];
   end;

   % find the common brain voxels that exist for all the runs
   %
   coord_mat = [];
   for i=run_idx,

      % load the coords and dims for each run
      %
      run_info = session_info.run(i);
      datamat_prefix = session_info.datamat_prefix;
      datamat_file = sprintf('%s_run%d.mat',datamat_prefix,i);
      out_name = fullfile(pls_data_path,datamat_file);
      eval(['load ' out_name ' coords dims']);

      c = zeros(1,prod(dims));
      c(coords) = 1;

      coord_mat = [coord_mat; c];
   end;

   coords_info.common = find(sum(coord_mat,1) == length(run_idx));
   coords_info.dims = dims;

   % compute the indices of the individual coords to get the common coords
   %
   coords_info.coords_idx = cell(1,length(run_idx));
   for i=1:length(run_idx),
      c = coord_mat(i,:); 
      run_coords = find(c == 1);
      c(coords_info.common) = c(coords_info.common) + 1;

      coords_info.coords_idx{i} = find(c(run_coords) == 2);
   end;
   

   if (nargout >= 1),
     coords_struct = coords_info;
   else
     try
        coords_file = sprintf('%s_common_coords.mat',datamat_prefix);
        coords_file = fullfile(pls_data_path,coords_file);
        save(coords_file,'coords_info','run_idx');
     catch
        errmsg = sprintf('ERROR: Cannot save coords info to the file: \n   %s.',coords_file);
        errordlg(errmsg,'Save Common Coords Error');
	waitfor(gcf);
        return;
     end;

     msg = ['** The common coords has been saved in ' coords_file];
     ShowProgress(progress_hdl,msg);
   end;

   return;


%-------------------------------------------------------------------------
function hdl = ShowProgress(progress_hdl,info)

  %  'initialize' - return progress handle if any
  %
  if ischar(progress_hdl) & strcmp(lower(progress_hdl),'initialize');
     if ~isempty(gcf) & isequal(get(gcf,'Tag'),'ProgressFigure'),
         hdl = gcf;
     else
         hdl = [];
     end;
     return;
  end;

  if ~isempty(progress_hdl)
     if ischar(info)
         rri_progress_status(progress_hdl,'Show_message',info);
     else
         rri_progress_status(progress_hdl,'Update_bar',info);
     end;
     return;
  end;

  if ischar(info),
     disp(info)
  end;

  return;                                       % ShowProgress

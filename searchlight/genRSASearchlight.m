function [smm_rs, smm_ps, n, searchlightRDMs] = genRSASearchlight(fullBrainVolumes, models, mask, userOptions, localOptions)

	% ARGUMENTS
	% fullBrainVolumes	A voxel x condition x session matrix of activity
	% 				patterns.
	%
	% models		A struct of model RDMs.
	%
	% mask     		A 3d or 4d mask to perform the searchlight in.
	%
	% userOptions and localOptions
	%
	% RETURN VALUES
	% smm_rs        4D array of 3D maps (x by y by z by model index) of
	%               correlations between the searchlight pattern similarity
	%               matrix and each of the model similarity matrices.
	%
	% smm_ps        4D array of 3D maps (x by y by z by model index) of p
	%               values computed for each corresponding entry of smm_rs.
	%
	% n             an array of the same dimensions as the volume, which
	%               indicates for each position how many voxels contributed
	%               data to the corresponding values of the infomaps.
	%               this is the number of searchlight voxels, except at the
	%               fringes, where the searchlight may illuminate voxels
	%               outside the input-data mask or voxel with all-zero
	%               time-courses (as can arise from head-motion correction).
	%
	% mappingMask_actual
	%               3D mask indicating locations for which valid searchlight
	%               statistics have been computed.
	%
	% Based on Niko Kriegeskorte's searchlightMapping_RDMs.m
	%
	% Additions by Cai Wingfield 2-2010:
	% 	- Now skips points in the searchlight where there's only one voxel inside.
	% 	- Now takes a userOptions struct for the input parameters.

	localOptions = setIfUnset(localOptions, 'averageSessions', true);

	%% Figure out whether to average over sessions or not
	if localOptions.averageSessions
		for sessionNumber = 1:size(fullBrainVolumes,3)
			thisSessionId = ['s' num2str(sessionNumber)];
			t_patsPerSession.(thisSessionId) = fullBrainVolumes(:,:,sessionNumber)';
		end%for:sessionNumber
	else
		justThisSession = 1;
		t_pats = fullBrainVolumes(:,:,justThisSession)';
		
		fprintf(['\nYou have selected not to average over sessions.\n         Only session number ' num2str(justThisSession) ' will be used.\n']);
		
	end%if

	%% Get parameters
	voxSize_mm = userOptions.voxelSize;
	searchlightRad_mm = userOptions.searchlightRadius;
	monitor = localOptions.monitor;
	nConditions = size(fullBrainVolumes, 2);
	
	clear fullBrainVolumes;

	% Prepare models
	modelRDMs_ltv = permute(unwrapRDMs(vectorizeRDMs(models)), [3 2 1]);

	% Prepare masks
	mask(isnan(mask)) = 0; % Just in case!
	if ndims(mask)==3
		inputDataMask=logical(mask);
		mappingMask_request=logical(mask);
	else
		inputDataMask=logical(mask(:,:,:,1));
		mappingMask_request=logical(mask(:,:,:,2));
	end

	% Check to see if there's more data than mask...
	if localOptions.averageSessions
		for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
			thisSessionId = ['s' num2str(sessionNumber)];
			t_patsPerSession.(thisSessionId) = t_patsPerSession.(thisSessionId)(:, inputDataMask(:));
		end%for:sessionNumber
	else
		if (size(t_pats,2)>sum(inputDataMask(:)))
			t_pats=t_pats(:,inputDataMask(:));
		end%if
	end%if

	% Other data
	volSize_vox=size(inputDataMask);
	nModelRDMs=size(modelRDMs_ltv,1);
	rad_vox=searchlightRad_mm./voxSize_mm;
	minMargin_vox=floor(rad_vox);


	%% create spherical multivariate searchlight

    [ctrRelSphereSUBs, nSearchlightVox] = make_searchlight_sphere(searchlightRad_mm, voxSize_mm, monitor);

	%% define masks
    [mappingMask_request_INDs, nVox_mappingMask_request, volIND2YspaceIND] = make_valid_mask(inputDataMask, localOptions);

	%% Initialize output
    %  n voxels contributing to infobased t at each location
	n=nan(volSize_vox);

	% similarity-graph-map the volume with the searchlight
	smm_bestModel=nan(volSize_vox);
	smm_ps=nan([volSize_vox,nModelRDMs]);
	smm_rs=nan([volSize_vox,nModelRDMs]);
	searchlightRDMs = nan([nConditions, nConditions, volSize_vox]);

	if monitor
		h_progressMonitor=progressMonitor(1, nVox_mappingMask_request,  'Similarity-graph-mapping...');
	end

	%% THE BIG LOOP! %%

	for cMappingVoxI=1:nVox_mappingMask_request
		
		if mod(cMappingVoxI,1000)==0
			if monitor
				progressMonitor(cMappingVoxI, nVox_mappingMask_request, 'Searchlight mapping Mahalanobis distance...', h_progressMonitor);
				%                 cMappingVoxI/nVox_mappingMask_request
			else
				fprintf('.');
			end%if
		end%if

		[x y z]=ind2sub(volSize_vox,mappingMask_request_INDs(cMappingVoxI));

		% compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;

		% exclude out-of-volume voxels
		outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
					cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
					cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));

		cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);

		% list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));

		% restrict searchlight to voxels inside validDataBrainMask
		cIllValidVox_volINDs=cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
		cIllValidVox_YspaceINDs=volIND2YspaceIND(cIllValidVox_volINDs);

		% note how many voxels contributed to this locally multivariate stat
		n(x,y,z)=length(cIllValidVox_YspaceINDs);
		
		if n(x,y,z) < 2, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
		
		if localOptions.averageSessions
			searchlightRDM = zeros(localOptions.nConditions, localOptions.nConditions);
			for session = 1:localOptions.nSessions
				sessionId = ['s' num2str(session)];
				searchlightRDM = searchlightRDM + squareform(pdist(t_patsPerSession.(sessionId)(:,cIllValidVox_YspaceINDs),'correlation'));
			end%for:sessions
			searchlightRDM = searchlightRDM / localOptions.nSessions;
		else
			searchlightRDM = squareform(pdist(t_pats(:,cIllValidVox_YspaceINDs), 'correlation'));
		end%if
		
		searchlightRDM = vectorizeRDM(searchlightRDM);
		
		% Locally store the full brain's worth of indexed RDMs.
		searchlightRDMs(:, :, x, y, z) = squareform(searchlightRDM);
		
		try
			[rs, ps] = corr(searchlightRDM', modelRDMs_ltv', 'type', 'Spearman', 'rows', 'pairwise');
		catch
			[rs, ps] = corr(searchlightRDM', modelRDMs_ltv, 'type', 'Spearman', 'rows', 'pairwise');
		end%try
		
		if localOptions.fisher
			for i = 1:numel(rs)
				rs(i) = fisherTransform(rs(i));
			end%for:i
		end%if
		
	%	[ignore, bestModelI] = max(rs);
		
	%    smm_bestModel(x,y,z) = bestModelI;
		smm_ps(x,y,z,:) = ps;
		smm_rs(x,y,z,:) = rs;
		
	end%for:cMappingVoxI

	%% END OF THE BIG LOOP! %%

	if monitor
		fprintf('\n');
		close(h_progressMonitor);
	end

	mappingMask_actual=mappingMask_request;
	mappingMask_actual(isnan(sum(smm_rs,4)))=0;

	%% visualize
	if monitor
		aprox_p_uncorr=0.001;
		singleModel_p_crit=aprox_p_uncorr/nModelRDMs; % conservative assumption model proximities nonoverlapping

		smm_min_p=min(smm_ps,[],4);
		smm_significant=smm_min_p<singleModel_p_crit;

		vol=map2vol(mask);
		vol2=map2vol(mask);
		
		colors=[1 0 0
				0 1 0
				0 1 1
				1 1 0
				1 0 1];
		
		for modelRDMI=1:nModelRDMs
			vol=addBinaryMapToVol(vol, smm_significant&(smm_bestModel==modelRDMI), colors(modelRDMI,:));
	% 		vol2=addBinaryMapToVol(vol2, smm_bestModel==modelRDMI, colors(modelRDMI,:));
		end
		
		showVol(vol);
		
	% 	vol2 = vol2*0.1;
	% 	vol2(vol<1) = vol(vol<1);
	% 	
	% 	showVol(vol2);
		
	end%if

end%function
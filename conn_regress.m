function conn_regress(seed_timeseries, motion_params, regressors, output_base)
% conn_regress.m 
%		This script regresses out the selected effects from the seed region 
%		and outputs text files that can be used in later voxel-wise regression
%		models.
%
% Usage: conn_regress(seed_timeseries, motion_params, regressors, output_base)
%
%		-- required arguments --
%		seed_timeseries     file that contains the vector for the seed region
%		motion_params       file that contains the motion parameters
%		regressors          cell array of filenames that contain additional regressors
%                                   required argument, but can be an empty cell array '{}' if desired.
%		output_base         a base string from which output filenames will be built
%
% Outputs: 
%		covariate file      Contains the covariates that were regressed from the
%                       seed region. The file is named for the subject and 
%                       the covariates included.
%		residuals file      Contains the residuals from the seed region after the
%                       covariates and linear trend have been regressed out.
%                       The file is named for the subject and _residuals_ and
%                       the covariates that were regressed out.
%

%  There is an internal check to make sure that the input files actually
%  exist. If the files don't exist, you will see the following message:
%  
%             Argument 'xxxxxx' failed validation
%
% Version 1.0 10/8/2009
% Created by Donald McLaren (dm10@medicine.wisc.edu) and Kim Farbota
% (farbota@wisc.edu).
% Wisconsin Alzheimer's Disease Center Imaging Core (brainmap.wisc.edu)
%
	%% Input Parsing
	p = inputParser;
	p.addRequired('seed_timeseries', @(x)exist(x,'file') && ~isempty(strfind(x,'.txt')));
	p.addRequired('motion_params', @(x)exist(x,'file'));
	p.addRequired('regressors', @(x)files_exist(x));
	p.addRequired('output_base', @(x)~exist(x,'file'));
	try
		p.parse(seed_timeseries, motion_params, regressors, output_base);
	catch exception
		disp(['ERROR parsing: ' exception.message]); disp(' ');
		help conn_regress
		return
	end
	P=p.Results;

	%% load the timeseries data from the files
	voi = multi_load(P.seed_timeseries);
	motion = multi_load(P.motion_params);
	regressors = multi_load(P.regressors);
	
	%% Ensure that all our timeseries have the same number of timepoints
	series_length = length(voi.data);
	validate_series_length(motion, series_length);
	for i = 1:length(regressors)
		validate_series_length(regressors(i), series_length);
	end
	
	%% Compute derivatives and compile covariates
  covariates = [motion.data fodiff(motion.data)];
	for i = 1:length(regressors)
		covariates = [covariates regressors(i).data fodiff(regressors(i).data)];
	end

  %% Create linear terms for regression
  trend(:,1) = transpose(-1:2/(series_length-1):1);
  trend(:,2) = 1;

  %% Compute residuals
  [B, BINT, R] = regress(voi.data, [covariates trend]);

	%% Write results to disk
	cov_tag = 'covar';
	resid_tag = 'resid';
	cov_file = [P.output_base '_' cov_tag '.txt'];
	resid_file = [P.output_base '_' resid_tag '.txt'];
	dlmwrite(cov_file, covariates, ' ');
  dlmwrite(resid_file, R, ' ');
end




%%%% AUXILIARY FUNCTIONS %%%%

function difference=fodiff(data)
% Computes first order difference on a matrix where rows are the timepoints
	difference = data-[0*data(end,:); data(1:end-1,:)];
end

function validate_series_length(series,proper_length)
% checks that number of timepoints the the series matches the given paramter
	length_err = 'ERROR: Incorrect series length in ';
	if (length(series.data) ~= proper_length)
		error([length_err series.filename]);
	end
end

function answer=files_exist(filenames)
% true iff all filenames in the given cell array exist.
	answer = all(cellfun(@exist,filenames));
end

function loads=multi_load(filenames)
% loads all files in the given cell array and returns a cell array of the 
% loaded contents
	loads = struct([]);
	if ischar(filenames), filenames = { filenames };, end;
	
	if iscellstr(filenames)
		for i = 1:length(filenames)
			loads(i).filename = filenames{i};
			loads(i).data = load(filenames{i});
		end
	else
		error('ERROR: multi_load(filenames) takes either a string or a cell array of strings');
	end
end

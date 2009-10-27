function Connectivity_regress(subject,region,timeseries,motion,varargin)
% Connectivity_regress.m 
%       This script regresses out the selected effects from the seed region and
%       outputs text files that can be used in later voxel-wise regression
%       models.
%
% Usage: iptest(subject, region, timeseries, motion, [options])
%
%       -- required arguments --
%       subject             unique subject identifier
%       region              name of the seed region for the analysis
%       timeseries          file that contains the vector for the seed region
%       motion              file that contains the motion parameters
%
%       -- optional regressors --
%       whitematter, FILE   white matter signal
%       ventricles, FILE    ventricle signal
%       brain, FILE         global brain signal
%       task, FILE          task-related regressor
%       cardiac, FILE       cardiac or other physiological noise
%       partials, FILE      signal from other seed regions to partial out
%
% Example: R=Connectivity_regress('2041','DMPFC','DMPFC_ts.txt','rp1.txt.', 'cardiac','cardiac.txt')
%
% Outputs: 
%       covariate file      Contains the covariates that were regressed from the
%                           seed region. The file is named for the subject and 
%                           the covariates included.
%       residuals file      Contains the residuals from the seed region after the
%                           covariates and linear trend have been regressed out.
%                           The file is named for the subject and _residuals_ and
%                           the covariates that were regressed out.
%

%  There is an internal check to make sure that the input files actually
%  exist. If the files don't exist, you will see the following message:
%  
%             Argument 'xxxxxx' failed validation
%
% Version 1.0 10/8/2009
% Created by Donald McLaren (dm10@medicine.wisc.edu) and Kim Farbota
% (farbota@wisc.edu).
% Wisconsin Alzheimer's Disease Center Imaging Core (www.brainmap.wisc.edu)
%

    %% Input Parsing
    p = inputParser;
    p.addRequired('subject', @ischar);
    p.addRequired('region', @ischar);
    p.addRequired('timeseries', @(x)exist(x,'file') && ~isempty(strfind(x,'.txt')));
    p.addRequired('motion', @(x)exist(x,'file'));
    p.addOptional('whitematter', [], @(x)exist(x,'file'));
    p.addOptional('ventricles', [], @(x)exist(x,'file'));
    p.addOptional('brain',[], @(x)exist(x,'file'));
    p.addOptional('task',[], @ischar);
    p.addOptional('cardiac',[], @(x)exist(x,'file'));
    p.addOptional('partials',[], @(x)exist(x,'file'));
    try
        p.parse(subject, region, timeseries, motion, varargin{:});
    catch exception
        disp(['ERROR parsing: ' exception.message]); disp(' ');
        help Connectivity_regress
        return
    end
    P=p.Results;

    %% Can do either cardiac or ventricles/whitematter but not both
    if ~xor(any(P.cardiac), (any(P.ventricles) && any(P.whitematter)))
        error('ERROR: cardiac OR ventricals/whitematter can be included, but not both or neither.')
    end

    %% Initialize Timeseries and Covariates
    VOI = load(P.timeseries);
    zero_vector = zeros(length(VOI),1);
    motion = zero_vector; wm = zero_vector; vent = zero_vector;
    brain = zero_vector; task = zero_vector; cardiac = zero_vector; partial = zero_vector;
        
    try
        motion = load(P.motion);
        if P.whitematter, wm = load(P.whitematter);, end;
        if P.ventricles, vent = load(P.ventricles);, end;
        if P.brain, brain = load(P.brain);, end;
        if P.task, task = load(P.task);, end;
        if P.cardiac, cardiac = load(P.cardiac);, end;
        if P.partials, partial = load(P.partials);, end;
    catch exception
        disp(['ERROR loading: ' exception.message]); disp(' ');
        return
    end

    if (P.partials)
        partialname = regexp(P.partials,'/(\w*)_rest','tokens');
        partialname = partialname((1));
    else
        partialname = '';
    end
    

    %% Length Error Checking
    molen = length(motion);
    length_err = 'Motion file length does not match the length of ';
    if (molen ~= length(vent)), error([length_err 'ventricle signal']), end;
    if (molen ~= length(wm)), error([length_err 'white matter signal']), end;
    if (molen ~= length(brain)), error([length_err 'global brain signal']), end;
    if (molen ~= length(task)), error([length_err 'task regressors']), end;
    if (molen ~= length(cardiac)), error([length_err 'cardiac regressors']), end;
    if (molen ~= length(partial)), error([length_err 'partial regressors']), end;

    %% Compute derivatives
    motion = [motion deriv(motion)];
    wm = [wm deriv(wm)];
    vent = [vent deriv(vent)];
    brain = [brain deriv(brain)];

    %% Create linear terms for regression
    trend(:,1) = transpose(-1:2/(length(motion)-1):1);
    trend(:,2) = 1;

    %% Compile all regressors together
    covariates = [motion wm vent brain task partial cardiac];
    covariates = covariates(:,any(covariates));

    %% Compute residuals
    [B, BINT, R] = regress(VOI, [covariates trend]);

    %% Output file tags
    motion_derivs_tag = 'motion_deriv';
    residuals_tag = 'residuals';
    cardiac_tag = 'cardiac';
    non_cardiac_tag = 'wm_vent';
    task_tag = 'task';
    brain_tag = 'brain';
    partial_tag = partialname;

    %% Assemble relevent output tags
    outtags = {};
    outtags(end+1) = { ternary(P.cardiac, cardiac_tag, non_cardiac_tag) };
    if P.task, outtags(end+1) = {task_tag}, end;
    if P.brain, outtags(end+1) = {brain_tag}, end;
    if P.partials, outtags(end+1) = {partials_tag}, end;

    %% Construct output file names
    covout = [subject '_' join(outtags,'_') '_' timeseries];
    residout = [subject '_' residuals_tag '_' join(outtags,'_') '_' timeseries];
    
    covout
    covariates
    
    residout
    R

    % dlmwrite(covout, covariates, ' ');
    % dlmwrite(residout, R, ' ');

end


function dmat = deriv(omat)
	dmat = omat-[0*omat(end,:); omat(1:end-1,:)];
end
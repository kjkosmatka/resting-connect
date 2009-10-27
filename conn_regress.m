function conn_regress(seed_timeseries, motion_params, regressors)
% Connectivity_regress.m 
%       This script regresses out the selected effects from the seed region and
%       outputs text files that can be used in later voxel-wise regression
%       models.
%
% Usage: conn_regress(seed_timeseries, motion_params, regressors={})
%
%       -- required arguments --
%       timeseries          file that contains the vector for the seed region
%       motion              file that contains the motion parameters
%       regressors          a cell array of filenames that contain additional regressors
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
    p.addRequired('seed_timeseries', @(x)exist(x,'file') && ~isempty(strfind(x,'.txt')));
    p.addRequired('motion_params', @(x)exist(x,'file'));
    p.addOptional('regressors',{}, @(x)files_exist(x));
    try
        p.parse(seed_timeseries, motion_params, regressors);
    catch exception
        disp(['ERROR parsing: ' exception.message]); disp(' ');
        help conn_regress
        return
    end
    P=p.Results;
    
    voi = load(P.seed_timeseries);
    motion = load(P.motion_params);
    regs = {}
    for i = 1:length(P.regressors)
        
    end
end


function result=files_exist(filenames)
result = all(cellfun(@exist,filenames));
end
function [tmsd,msd,msd_N,msd_std] = compute_MSD(x, dt, varargin)
% [tmsd,msd,msd_N,msd_std] = compute_MSD(x, dt, varargin)
%
% Compute the mean squared displacements of trajectories.
%
% Inputs:
% x: Trajectories sampled uniformly in time. Ncells by Ntime matrix.
% dt: Time step. scalar
%
% Optional inputs:
% max_tau: Maximum time delay over which to compute the MSD. scalar. default is the whole trajectory.
% every_n: Compute the MSD at every every_n time points. i.e. subsample the MSD. scalar. default value is 1, which computes the MSD at every time point.
%
% Outputs:
% tmsd: Time points at which the MSD was computed. 1 by N vector.
% msd: MSD of each trajectory. Ncells by N matrix.
% msd_N: Number of observations for each value in the matrix msd. Ncells by N matrix.
% msd_std: Standard deviation of the squared displacement samples. Ncells by N matrix.
%
% Henry H. Mattingly, November 2023

max_tau=(size(x,2)-1)*dt;
every_n=1;

if nargin>=3
    if ~isempty(varargin{1})
        max_tau=varargin{1};
    end
end
if nargin>=4
    if ~isempty(varargin{2})
        every_n=varargin{2};
    end
end
% etc

Ncells = size(x,1);
nt = size(x,3);

max_frames = 1+floor(max_tau/dt);
nt = min(nt, max_frames);

eval_inds = 1:every_n:nt;
neval = length(eval_inds);
tmsd = (eval_inds-1)*dt;

% compute vacf, etc for time delay less than or equal to max_tau
msd = nan(Ncells,neval);
msd_N = nan(Ncells,neval);
msd_std = nan(Ncells,neval);

for j = 1:length(eval_inds)
    tempx = sum((x(:,:,1:end-eval_inds(j)+1)-x(:,:,eval_inds(j):end)).^2,2);
    msd(:,j) = nanmean(tempx,3);
    msd_N(:,j) = sum(~isnan(tempx),3);
    msd_std(:,j) = nanstd(tempx,[],3);
    
end

% msd_se = msd_std;


% disp('done.')

end


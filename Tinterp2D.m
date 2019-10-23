% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/1/2019

%% Functionality
% This function calculates the mean vegetation index of a time step by averaging
%  its two nearest neighbors in time. The weighting factors are defined by the
%  time distances from the two neighbors to the time step.

%% Input
% vi: space-time class (V2DTCls.m) object for vegetation index before and after
%      the time step of interest;
% dn: date number of the time step of interest.

%% Output
% vi: vegetation index of the time step of interest.

%% Additional note
% Require V2DCls.m.

function VI=Tinterp2D(vi,dn)
%% Check inputs
narginchk(2,2);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'vi',@(x) validateattributes(x,{'double','V2DTCls'},{'nonempty'},mfilename,'vi'));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dn'));

parse(ips,vi,dn);
clear ips

%% Mean VI
[X,~]=vi.GridCls;

VI=[];
for n=1:length(vi.Fnm)
  VI=cat(2,VI,reshape(vi.readCls(n),numel(X),1));
end
VI(isnan(VI(:,1)),1)=VI(isnan(VI(:,1)),2); % Fill NaN by each other's records
VI(isnan(VI(:,2)),2)=VI(isnan(VI(:,2)),1);

Dn=vi.TimeCls;
w=1-abs(dn-Dn)/diff(Dn);
VI=reshape(VI*w,size(X));
end

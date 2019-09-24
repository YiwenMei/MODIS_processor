% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/01/2019

%% Functionality
% This function calculates the mean vegetation index of a time step by averaging
%  its two nearest neighbors in time. The weighting factors are defined by the
%  time distances from the two neighbors to the time step.
%% Input
% vib/: details of file or workspace variable for vegetation index before/after
%  via  the time step of interest;
%  dn : date number of the time step of interest.

%% Output
% vi: vegetation index of the time step of interest.

%% Additional note
% Require read2Dvar.m.

function vi=Tinterp2D(vib,via,dn)
%% Check inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 3 required inputs\n',mfilename);

addRequired(ips,'vib',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'vib',1));
addRequired(ips,'via',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'via',2));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dn',3));

%% Mean VI
[~,dnb,~]=fileparts(vib{1});
dnb=datenum(dnb(3:end),'yyyymmdd'); % datenum of VI record before the inqury time step
vib=read2Dvar(vib);
wb=1-abs((dn-dnb))/8;

[~,dna,~]=fileparts(via{1});
dna=datenum(dna(3:end),'yyyymmdd'); % datenum of VI record after the inqury time step
via=read2Dvar(via);
wa=1-abs((dn-dna))/8;

vi=nansum(cat(3,vib*wb,via*wa),3); % mean of VI

k1=isnan(vib) | isnan(via);
vi(k1)=nansum([vib(k1) via(k1)],2);
vi(isnan(vib) & isnan(via))=NaN;
end

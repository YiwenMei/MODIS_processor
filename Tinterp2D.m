% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/01/2019

%% Functionality
% This function calculates the values of a variable for a time step by linearly
%  interpolating its two nearest neighbors in time.

%% Input
% imb/: details of file or matlab workspace variable for the variable before/after
%  ima  the time step of interest;
%  dn : date number of the time step of interest.

%% Output
% im_dn: variable values for the time step of interest.

%% Additional note
% Require read2Dvar.m.

function im_dn=Tinterp2D(imb,ima,dn)
%% Check inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 3 required inputs\n',mfilename);

addRequired(ips,'vib',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'vib',1));
addRequired(ips,'via',@(x) validateattributes(x,{'double','cell'},{'nonempty'},mfilename,'via',2));
addRequired(ips,'dn',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'dn',3));

%% Mean value
[~,dnb,~]=fileparts(imb{1});
dnb=datenum(dnb(3:end),'yyyymmdd'); % datenum of record before the inqury time step
imb=read2Dvar(imb);
wb=1-abs((dn-dnb))/8;

[~,dna,~]=fileparts(ima{1});
dna=datenum(dna(3:end),'yyyymmdd'); % datenum of record after the inqury time step
ima=read2Dvar(ima);
wa=1-abs((dn-dna))/8;

im_dn=nansum(cat(3,imb*wb,ima*wa),3);

k1=isnan(imb) | isnan(ima);
im_dn(k1)=nansum([imb(k1) ima(k1)],2);
im_dn(isnan(imb) & isnan(ima))=NaN;
end

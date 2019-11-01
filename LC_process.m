% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/7/2019

%% Functionality
% This function contains several functionalities:
%  1)read all MCD12Q1 LC Type 3 images for the study area (MODISimg.m);
%  2)mosaic all LC images together (MODISimg.m);
%  3)project and crop the mosaiced LC image (MODISimg.m);
%  4)calculate the fraction of a land cover class for a block (LC_Frac.m);
%  5)output the land cover fraction images (matV2tif,m).

%% Input
% lcfl : full name list of the land cover files;
%  tid : ID of the land cover type (possible value is 1 to 5, Friedl et al. 2010);
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: output directory for the images (e.g. C:\...\LCFrac\);
% Bound: boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;
%  ors : output coordinate system (e.g. 'EPSG:102012');

% thr : a percentage of nodata value for a block where block with this percent
%       of nodata value is assigned the nodata value;
% ndv : no-data value assigned to the output images.
% pflg: parallel flag (true - default, use parallel channel; false - squential);

%% Output
% LC: mosaiced land cover image of the study area at original resolution;

% Land cover fraction images are stored in oupth with the naming convention as
%  LCXX_yyyy.tif where XX represents the ID of a land cover class.

%% Additional Note
% 1)Index for the land cover classes of different types (Friedl et al. 2010):
%   Type 1 (IGBP)
%    01 - Evergreen Needleleaf     02 - Evergreen Broadleaf 03 - Deciduous Needleleaf
%    04 - Deciduous Broadleaf      05 - Mixed Forests       06 - Closed Shrublands
%    07 - Open Shrublands          08 - Woody Savannas      09 - Savannas
%    10 - Grasslands               11 - Permanent Wetlands  12 - Croplands
%    13 - Urban and Built-up Lands 14 - Cropland/Natural Vegetation Mosaics
%    15 - Permanent Snow and Ice   16 - Barren              17 - Water Bodies
%   Type 2 (UMD)
%    00 - Water Bodies         01 - Evergreen Needleleaf 02 - Evergreen Broadleaf
%    03 - Deciduous Needleleaf 04 - Deciduous Broadleaf  05 - Mixed Forests
%    06 - Closed Shrublands    07 - Open Shrublands      08 - Woody Savannas
%    09 - Savannas             10 - Grasslands           11 - Permanent Wetlands
%    12 - Croplands            13 - Urban and Built-up Lands
%    14 - Cropland/Natural Vegetation Mosaics            15 - Barren
%   Type 3 (LAI/FPAR)
%    00 - Water Bodies        01 - Grassland            02 - Shrublands
%    03 - Broadleaf Croplands 04 - Savannas             05 - Evergreen Broadleaf
%    06 - Deciduous Broadleaf 07 - Evergreen Needleleaf 08 - Deciduous Needleleaf
%    09 - Unvegetated         10 - Urban and Built-up Lands
%   Type 4 (BGC)
%    00 - Water Bodies         01 - Evergreen Needleleaf 02 - Evergreen Broadleaf
%    03 - Deciduous Needleleaf 04 - Deciduous Broadleaf  05 - Annual Broadleaf
%    06 - Annual Grass         07 - Non-Vegetated land   08 - Urban and Built-up Lands
%   Type 5 (PFT)
%    00 - Water Bodies             01 - Evergreen Needleleaf   02 - Evergreen Broadleaf
%    03 - Deciduous Needleleaf     04 - Deciduous Broadleaf    05 - Shrub
%    06 - Grass                    07 - Cereal Crop            08 - Broadleaf Crop
%    09 - Urban and Built-up Lands 10 - Permanent Snow and Ice 11 - Barren
% 2)Require MODISimg.m, LCFrac.m, and matV2tif.m;
% 3)Need to have GDAL installed to run the code.

function LC=LC_process(lcfl,tid,wkpth,oupth,Bound,ors,varargin)
%% Check the input
narginchk(6,9);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'lcfl',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'lcfl'));
addRequired(ips,'tid',@(x) validateattributes(x,{'double'},{'scalar','>=',1,'<=',5},mfilename,'tid'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'oupth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'oupth'));
addRequired(ips,'Bound',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Bound'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));

addOptional(ips,'thr',.5,@(x) validateattributes(x,{'double'},{'scalar','>=',0,'<=',1},mfilename,'thr'));
addOptional(ips,'ndv',-999,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,lcfl,tid,wkpth,oupth,Bound,ors,varargin{:});
thr=ips.Results.thr;
ndv=ips.Results.ndv;
pflg=ips.Results.pflg;
clear ips

%% Properties of input records
hif=hdfinfo(lcfl(1,:));
rn=hif.Vgroup.Name; % Name of the record
hif=hif.Vgroup.Vgroup(1).SDS(tid); % LC type
lcc=hif.Attributes(2:length(hif.Attributes)-3); % Land cover classes
vrf=hif.Name; % Name of the field
ndv_o=double(hif.Attributes(length(hif.Attributes)-2).Value); % no-data-value of LC
% lcc=hif.Attributes(3:length(hif.Attributes)-2); % Land cover classes
% ndv_o=double(hif.Attributes(length(hif.Attributes)-1).Value); % no-data-value of LC

%% Image processing (read, mosaic, project, crop, resample)
LC=MODISimg(lcfl,wkpth,'',[Bound(1:2,:);[NaN NaN]],ors,rn,vrf,'near',pflg);
LC(LC==ndv_o)=NaN;

%% Calculate land cover fraction
[~,ys,~]=fileparts(lcfl(1,:));
ys=cell2mat(regexp(ys,'.A(\d{4})\d{3}.','tokens','once'));

rsx_in=(Bound(2,1)-Bound(1,1))/size(LC,2);
rsy_in=(Bound(1,2)-Bound(2,2))/size(LC,1);

% The 1st tile for idn
idn=fullfile(wkpth,sprintf('idLC_r%i.mat',Bound(3,1)));
[~]=LC_process_sub(LC,ndv,Bound,idn,lcc,1,thr,rsx_in,rsy_in,oupth,ys,ors,wkpth);

% Process the other tiles
switch pflg
  case true
    parfor v=2:length(lcc)
      [~]=LC_process_sub(LC,ndv,Bound,idn,lcc,v,thr,rsx_in,rsy_in,oupth,ys,ors,wkpth);
    end

  case false
    for v=2:length(lcc)
      [~]=LC_process_sub(LC,ndv,Bound,idn,lcc,v,thr,rsx_in,rsy_in,oupth,ys,ors,wkpth);
    end
end
end

function LCF=LC_process_sub(LC,ndv,Bound,idn,lcc,v,thr,rsx_in,rsy_in,oupth,ys,ors,wkpth)
% Calculate the land cover fraction
LCF=LC_Frac(LC,ndv,Bound(:,1),Bound(:,2),idn,double(lcc(v).Value),thr,...
    [Bound(1:2,1);rsx_in],[Bound(1:2,2);rsy_in]);

% Output the land cover fraction images
IMo=fullfile(oupth,sprintf('LC%02i_%s.tif',lcc(v).Value,ys));
matV2tif(IMo,LCF,Bound(1,1),Bound(2,2),Bound(3,1),ndv,ors,wkpth);
end

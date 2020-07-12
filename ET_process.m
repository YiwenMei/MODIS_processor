% Yiwen Mei (ymei@gmu.edu)
% CEIE, George Mason University
% Last update: 4/28/2020

%% Functionality
% This function includes several functionalities:
%  1)read and average the variable of interest (ET/PET/LE/PLE) from MOD16A2 and
%     MYD16A2 for each tile;
%  2)convert ET/PET to mm/day and output the averaged variable as geotiff (doy2date.m
%     and matV2tif.m);
%  3)mosaic all the tiles files together (MODISimg.m);
%  4)project, crop and resample the mosaiced variable (MODISimg.m).

%% Input
% ETfl : a 2-by-1 cell array storing the name lists of .hdf files for MOD16A1
%         and MYD11A1 (in cell 1 and 2, repectively) for all the tiles of a time
%         step;
%  vid : index for the variable within the .hdf file (1 - ET_500m, 2 - LE_500m,
%         3 - PET_500m, 4 - PLE_500m, & 5 - ET_QC_500m);
% wkpth: working directory for the code (e.g., /.../wkdir/);
% opth : output directory for the images (e.g., /.../ET/);
%  GIf : boundary latitude and longitude and resolution of output images (it
%         follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%         l/r/b/t stands for left/right/bottom/top;
%  ors : output coordinate system (e.g., 'EPSG:102012');
%  ndv : no-data value assigned to the output files.

% pflg: parallel flag (true - default, use parallel channel; false - squential);

%% Output
% Output images are stored in opth as MODIS-<variable name>-yyyymmdd.tif (variable
%  names can be ET, LE, PET, PLE, or ET_QC based on the vid used).

%% Additional note
% Require MODISimg.m, doy2date.m, and matV2tif.m.

function imo=ET_process(ETfl,vid,wkpth,opth,GIf,ors,ndv,varargin)
%% Check the input
narginchk(7,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'ETfl',@(x) validateattributes(x,{'cell'},{'vector'},mfilename,'ETfl'));
addRequired(ips,'vid',@(x) validateattributes(x,{'double'},{'scalar','>=',1,'<=',5},mfilename,'vid'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'GIf'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));
addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));

addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,ETfl,vid,wkpth,opth,GIf,ors,ndv,varargin{:});
pflg=ips.Results.pflg;
clear ips varargin

%% Properties of input records
ETfl=cell2mat(ETfl); % Check whether the list is empty
% ors_o='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"';
ors_o='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"';

if ~isempty(ETfl)
  hif=hdfinfo(ETfl(1,:));
  vrf=hif.Vgroup.Vgroup(1).SDS(vid).Name;
  hif=hif.Vgroup.Vgroup(1).SDS(vid);
  scf=double(hif.Attributes(5).Value); % scale factor
  ndv_o=double(hif.Attributes(3).Value); % no-data value of ET

%% Find unique tiles
  A=cell(size(ETfl,1),1);
  for t=1:size(ETfl,1)
    [a,~]=regexp(ETfl(t,:),'h+\d{2}v+\d{2}','match');
    A{t}=cell2mat(a);
  end
  A=unique(A);

%% Processing (read, mosaic, project, crop, resample)
  Tfl=[];
  ETfl=mat2cell(ETfl,ones(size(ETfl,1),1),size(ETfl,2));
  switch pflg
    case true
      parfor n=1:length(A)
        tfn=ET_process_sub(ETfl,A{n},vrf,ndv_o,scf,ndv,ors_o,wkpth);
        Tfl=[Tfl;tfn];
      end

    case false
      for n=1:length(A)
        tfn=ET_process_sub(ETfl,A{n},vrf,ndv_o,scf,ndv,ors_o,wkpth);
        Tfl=[Tfl;tfn];
      end
  end
  tfl=mat2cell(Tfl,ones(size(Tfl,1),1),size(Tfl,2));
  tfl=regexp(tfl,'(?<year>\d{4})(?<day>\d{3})','match');
  ds=cellfun(@cell2mat,tfl,'UniformOutput',false);
  ds=cell2mat(unique(ds));
  ds=datestr(doy2date(str2double(ds(5:end)),str2double(ds(1:4))),'yyyymmdd');
  a=regexp(vrf,'_\d+m');
  oun=fullfile(opth,sprintf('MODIS-%s-%s.tif',vrf(1:a-1),ds));
  imo=MODISimg(Tfl,wkpth,oun,GIf,ors);
  imo(imo==ndv)=NaN;
end
end

function tfn=ET_process_sub(ETfl,a,vrf,ndv_o,scf,ndv,ors_o,wkpth)
hfl=ETfl(contains(ETfl,a));
[~,ds,~]=fileparts(hfl{1});
ds=cell2mat(regexp(ds,'A(\d{7}).h','tokens','once'));

% Number of accumulation day
switch vrf
  case {'ET_500m','PET_500m'} % kg/m^2/day = mm/day
    if ~isleap(str2double(ds(1:4))) && strcmp(ds(5:7),'361')
      cdl=5;
    elseif isleap(str2double(ds(1:4))) && strcmp(ds(5:7),'361')
      cdl=6;
    else
      cdl=8;
    end
  otherwise % J/m^2/day
    cdl=1;
end

% Average the MOD and MYD record
Et=[];
for t=1:length(hfl)
  et=double(hdfread(hfl{t},vrf));
  et(et>max(ndv_o) | et<min(ndv_o))=NaN;
  Et=cat(3,Et,et);
end
Et=nanmean(Et,3)*scf/cdl;
Et(isnan(Et))=ndv;

% Output the averaged ET of a tile
hif=hdfinfo(hfl{1},'eos');
xll=hif.Grid.UpperLeft(1,1);
yll=hif.Grid.LowerRight(1,2);
rs_o=mean(abs(hif.Grid.UpperLeft-hif.Grid.LowerRight)./[hif.Grid.Columns hif.Grid.Rows]);
tfn=fullfile(wkpth,sprintf('%s-%s-%s.tif',vrf,ds,a));
matV2tif(tfn,Et,xll,yll,rs_o,ndv,ors_o,wkpth);
end

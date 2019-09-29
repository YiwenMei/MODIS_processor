% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/9/2019

%% Functionality
% This function includes several functionalities:
%  1)read all MOD13Q1/MYD13Q1 NDVI/EVI tiles for the study area (resizeimg.m and matV2tif.m);
%  2)mosaic all NDVI/EVI tiles together (MODISimg.m);
%  3)project, crop, resample, and output the mosaiced VI image (MODISimg.m and doy2date.m);

%% Input
% vifl : full name list of the MODIS vegetation index files;
%  vrf : field name of variable (e.g. '250m 16 days NDVI');
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: path to store the output VI images (e.g. C:\...\VegIdx\);
% Bound: boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;
%  ors : output coordinate system (e.g. 'EPSG:102012');

% thr : a percentage of nodata value for a block where block with this percent
%       of nodata value is assigned the nodata value;
% ndv : no-data value assigned to the output images.
% pflg: parallel flag (true - default, use parallel channel; false - squential);

%% Output
% VI: processed vegetation index image read as Matlab variable of the study area;

% Output is the processed vegetation index image stored in oupth. The naming
%  convention is VIyyyymmdd.tif.

%% Additional Note
% Require MODISimg.m, doy2date.m, resizeimg.m, and matV2tif.m.

function VI=VI_process(vifl,vrf,wkpth,oupth,Bound,ors,varargin)
%% Check the input
narginchk(6,9);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 6 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'vifl',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vifl',1));
addRequired(ips,'vrf',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vrf',2));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth',3));
addRequired(ips,'oupth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'oupth',4));
addRequired(ips,'Bound',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Bound',5));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors',6));

addOptional(ips,'thr',.5,@(x) validateattributes(x,{'double'},{'scalar','>=',0,'<=',1},mfilename,'thr',7));
addOptional(ips,'ndv',-999,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv',8));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg',9));

parse(ips,vifl,vrf,wkpth,oupth,Bound,ors,varargin{:});
thr=ips.Results.thr;
ndv=ips.Results.ndv;
pflg=ips.Results.pflg;
clear ips

%% Properties of input records
hif=hdfinfo(vifl(1,:));
hif=hif.Vgroup.Vgroup(1).SDS(1);
scf=double(hif.Attributes(5).Value); % scale factor
ndv_o=double(hif.Attributes(4).Value); % no-data-value of LC

hif=hdfinfo(vifl(1,:),'EOS');
hif=hif.Grid;
xli=hif.UpperLeft(1);
xri=hif.LowerRight(1);
yti=hif.UpperLeft(2);
ybi=hif.LowerRight(2);
rxi=(xri-xli)/hif.Columns;
ryi=(yti-ybi)/hif.Rows;

kx=max(1,fix(Bound(3,1)/rxi));
ky=max(1,fix(Bound(3,2)/ryi));
if kx>1 || ky>1
%% Upscale the image
  Tfl=[];
  idn=fullfile(wkpth,['idVI_r' num2str(Bound(3,1),'%i') '.mat']);
  ors_i='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ';

% Write the upscaled image to the working directory
% The 1st tile for idn
  tfn=VI_process_sub(vifl(1,:),vrf,ndv_o,scf,ndv,kx,ky,idn,thr,wkpth,1,rxi,ors_i);
  Tfl=[Tfl;tfn];

% Process the other tiles
  switch pflg
    case true
      parfor n=2:size(vifl,1)
        tfn=VI_process_sub(vifl(n,:),vrf,ndv_o,scf,ndv,kx,ky,idn,thr,wkpth,n,rxi,ors_i);
        Tfl=[Tfl;tfn];
      end

    case false
      for n=2:size(vifl,1)
        tfn=VI_process_sub(vifl(n,:),vrf,ndv_o,scf,ndv,kx,ky,idn,thr,wkpth,n,rxi,ors_i);
        Tfl=[Tfl;tfn];
      end
  end

%% Process the upsacled image
  [~,ds,~]=fileparts(vifl(1,:));
  [ds,~]=regexp(ds,'.A(\d{4})(\d{3}).','tokens','once');
  ys=ds{1};
  ds=doy2date(str2double(ds{2}),str2double(ys));
  ds=datestr(ds,'yyyymmdd');

% Image processing (read, mosaic, project, crop, resample)
  oun=fullfile(oupth,['VI' ds '.tif']);
  VI=MODISimg(Tfl,wkpth,oun,Bound,ors,'','','bilinear',pflg);
end
end

function tfn=VI_process_sub(vfn,vrf,ndv_o,scf,ndv,kx,ky,idn,thr,wkpth,n,rxi,ors_i)
% Upscale vi
vi=double(hdfread(vfn,vrf));
vi(vi==ndv_o)=NaN;
vi=vi/scf;
if max(kx,ky)>1
  vi=resizeimg(vi,ndv,kx,ky,idn,thr);
end

% Write the upscaled image to the working directory
hif=hdfinfo(vfn,'EOS');
hif=hif.Grid;
xli=hif.UpperLeft(1);
ybi=hif.LowerRight(2);
[~,nm,~]=fileparts(vfn);
nm=cell2mat(regexp(nm,'A(?<year>\d+)(?<day>\d{3})','match'));
tfn=sprintf('%s%02i.tif',[wkpth nm '_p'],n);

matV2tif(tfn,vi,xli,ybi,kx*rxi,ndv,ors_i,wkpth);
end

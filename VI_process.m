% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/09/2018

%% Functionality
% This function includes several functionalities:
%  1)read all MOD13Q1/MYD13Q1 NDVI/EVI tiles for the study area;
%  2)upscale each tile to a coarser resolution (resizeimg.m);
%  3)mosaic all NDVI/EVI tiles together (MODISimg.m);
%  4)reproject, crop and resample the mosaiced VI image (MODISimg.m);

%% Input
% vifl : full name list of the MODIS vegetation index files;
%  vrf : field name of variable (e.g. '250m 16 days NDVI');
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: path to store the output VI images (e.g. C:\...\VegIdx\);
% rx/ry: resolution of output images;
% xl/xr: left/right side coordinate;
% yb/yt: same as xl/xr but for the bottom/top side's;
%  ors : output coordinate system (e.g. 'EPSG:102012');
%  thr : a percentage of nodata value for a block where block with this percent
%        of nodata value is assigned the nodata value;
%  ndv : no-data value assigned to the output images.

%% Output
% Output is a mosaiced vegetation index map stored in oupth. The naming convention
%  is VIyyyyddd.tif the same as the input date.

function VI=VI_process(vifl,vrf,wkpth,oupth,xl,xr,rx,yb,yt,ry,ors,thr,ndv)
%% Check the input
switch nargin
  case {1:12}; error('Not enough number of arguments');
  case 13
  otherwise; error('Too many number of arguments');
end

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

kx=max(1,fix(rx/rxi));
ky=max(1,fix(ry/ryi));
if kx>1 || ky>1
%% Upscale the image
  Tfl=[];
  idn=fullfile(wkpth,['idVI_r' num2str(rx,'%i') '.mat']);
  parfor n=1:size(vifl,1)
    vi=double(hdfread(vifl(n,:),vrf));
    vi(vi==ndv_o)=NaN;
    vi=vi/scf;
    vi=resizeimg(vi,ndv,kx,ky,idn,thr);

% Write the upscaled image to the working directory
    hif=hdfinfo(vifl(n,:),'EOS');
    hif=hif.Grid;
    xli=hif.UpperLeft(1);
    ybi=hif.LowerRight(2);

    vi_ors='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ';

    [~,nm,~]=fileparts(vifl(n,:));
    nm=regexp(nm,'A(?<year>\d+)(?<day>\d{3})','match');
    nm=cell2mat(nm);
    tfn=sprintf('%s%02i.tif',[wkpth nm '_p'],n);
    
    matV2tif(tfn,vi,xli,ybi,kx*rxi,ndv,vi_ors,wkpth);

    Tfl=[Tfl;tfn];
  end

%% Process the upsacled image
  [~,ns,~]=fileparts(vifl(1,:));
  ys=ns(10:13);
  ds=ns(14:16);
  ds=doy2date(str2double(ds),str2double(ys));
  ds=datestr(ds,'yyyymmdd');

% Image processing (read, mosaic, project, crop, resample)
  oun=fullfile(oupth,['VI' ds '.tif']);
  VI=MODISimg(Tfl,[],[],wkpth,oun,xl,xr,rx,yb,yt,ry,ors,'bilinear');
end
end

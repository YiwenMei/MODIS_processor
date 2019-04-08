% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 06/02/2018

%% Functionality
% This function includes several functionalities:
%  1)read all the tiles of MOD13Q1 vegetation index;
%  2)mask out the NDVI grids with snow/water covers (optional);
%  3)mosaic all the VI tiles;
%  4)reproject and crop the mosaiced VI.

%% Input
%   vfl  : strings of the list of MODIS vegetation index files
%   mfl  : strings of the list of MODIS mask files (set this to [] if functionality
%          # 2 is not required);
%  oupth : path to store the output NDVI geotiff
% rxo/ryo: resolution of interest in coor unit.
%  xl/xr : left/right corner coordinate in coor unit;
%  yb/yt : same as xl/xr but for the bottom and top corner coordinate;
%  coor  : coordinate system of interested;
%   ndv  : refer to resizeimg.m for more details;
%   thr  : refer to resizeimg.m for more details;

%% Output
% Output is a mosaiced vegetation index map stored in oupth. The naming convention
%  is VIyyyyddd.tif the same as the input date.

function mask_NDVI(vfl,mfl,oupth,rxo,ryo,xl,xr,yb,yt,coor,thr,ndv)
% Domain info
hi=nan(size(vfl,1),1);
vi=nan(size(vfl,1),1);
for n=1:size(vfl,1)
  vfn=vfl(n,:);
  [~,ns,~]=fileparts(vfn);
  hi(n)=str2double(ns(19:20));
  vi(n)=str2double(ns(22:23));
end
ntr=length(unique(vi)); % num. of tile in row
ntc=length(unique(hi)); % num. of tile in column
hi=min(hi)-1; % initial column tile num.
vi=min(vi)-1; % initial row tile num.

% VI tile info
hif=hdfinfo(vfn);
hif=hif.Vgroup.Vgroup(1).SDS(1);
ivn=hif.Name; % Name of VI
nr=hif.Dims(1).Size; % num. of grid in row of a tile
nc=hif.Dims(2).Size; % num. of grid in column of a tile
scf=double(hif.Attributes(5).Value); % scale factor of VI
ndv_vi=double(hif.Attributes(4).Value); % no-data value of VI
hif=hdfinfo(vfn,'eos');
hif=hif.Grid;
rtx=floor(rxo*hif.Columns/(hif.LowerRight(1)-hif.UpperLeft(1))); % Upscale ratio of x
rty=floor(ryo*hif.Rows/(hif.UpperLeft(2)-hif.LowerRight(2))); % Upscale ratio of y

% Mask tile info
if ~isempty(mfl)
  mfn=mfl(1,:);
  hif=hdfinfo(mfn);
  hif=hif.Vgroup.Vgroup(1).SDS(1); % 8-day max cover tile info
  iv13=hif.Name;
  nr1=hif.Dims(1).Size;
  nc1=hif.Dims(2).Size;

  mk_ndv=[0 1 11 50 254 255]; % Snow/water mask values
  mk_nv=25;
  mk_pv=[37 39 100 200];
end

%% Mosaic the VI tiles
VI=ndv*ones(ntr*nr/rty,ntc*nc/rtx);
xyll=nan(size(vfl,1),4);

for n=1:size(vfl,1)
% Read the VI
  vfn=vfl(n,:);
  ndvi=double(hdfread(vfn,ivn))/scf;

% Read the mask
  if ~isempty(mfl)
    mfn=mfl(n,:);

    mask=double(hdfread(mfn,iv13));
    k=multiFind(mask,mk_ndv);
    mask(k)=0;
    k=multiFind(mask,mk_nv);
    mask(k)=-1;
    k=multiFind(mask,mk_pv);
    mask(k)=1;
    mask=repelem(mask,nr/nr1,nc/nc1); % Nearest interpolate to the VI res

% Apply the mask
    ndvi(mask==1 | ndvi==ndv_vi/scf)=ndv;
    clear mask k
  end
  ndvi(ndvi==ndv_vi/scf)=ndv;

% Upscale VI
  hif=hdfinfo(vfn,'eos');
  hif=hif.Grid;
  xlt=hif.UpperLeft(1);
  xrt=hif.LowerRight(1);
  rx=(xrt-xlt)/hif.Columns; % resolution of a VI tile
  ybt=hif.LowerRight(2);
  ytt=hif.UpperLeft(2);
  ry=(ytt-ybt)/hif.Rows;
  ndvi=resizeimg(ndvi,xlt,xrt,rx,rtx,ybt,ytt,ry,rty,thr,ndv);

% Mosaic the NDVI tiles
  [~,ns,~]=fileparts(vfn);
  h=str2double(ns(19:20))-hi;
  v=str2double(ns(22:23))-vi;
  VI((v-1)*size(ndvi,1)+1:v*size(ndvi,2),(h-1)*size(ndvi,1)+1:h*size(ndvi,1))=ndvi;

  xyll(n,:)=[xlt ytt xrt ybt];
end
clear ndvi
%% Project the VI
% Write asc
ncn=[oupth 'p.asc'];

fid=fopen(ncn,'w');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(VI,2))],...
    ['nrows ' num2str(size(VI,1))],['xllcorner ' num2str(min(xyll(:,1)),12)],...
    ['yllcorner ' num2str(min(xyll(:,4)),12)],['cellsize ' num2str(rtx*rx,12)],...
    ['NODATA_value ' num2str(ndv)]);
dlmwrite(ncn,VI,'delimiter',' ','-append');
fclose(fid);

clear VI

% Project the VI
fun='gdalwarp -overwrite -of GTiff -r bilinear '; % GDAL function
pr1=['-t_srs ' coor ' '];
pr2=sprintf('%s %i %i %i %i ','-te',xl,yb,xr,yt);
pr3=sprintf('%s %i %i ','-tr',rxo,ryo);
pr4='-s_srs "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" ';
pr5=sprintf('%s %i ','-srcnodata',ndv);

par=[pr1 pr2 pr3 pr4 pr5];
ouv=[oupth 'VI' ns(10:16) '.tif'];
system([fun par '"' ncn '" "' ouv '"']);

delete(ncn);
end

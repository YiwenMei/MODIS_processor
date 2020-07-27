% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 4/20/2020

%% Functionality
% This function is used to process the GLEAM-daily product (Miralles et al. 2011).
%  Its main functionalities are
%   1)read one time step (one day) within the yearly file and
%   2)extract the variable record of the time step based on a given lat/lon box.
%  The extracted records will be outputted as .mat files to the output directory.
%  The optional functionalities of this function are
%   3)resample the extracted record;
%   4)project the extracted record to another coordinate system;
%   5)crop the projected record to a retangular box.
%  The operations are looped for all time steps contain in the yearly file and
%   the processed files are outputted as .tif in the output directory.

%% Input
% fname: full name with path of the input GLEAM daily record for a year (e.g.,
%         G:\GLEAM\E_2013_GLEAM_v3.3b.nc);
% wkpth: working directory of the code;
% opth : output directory to store the outputed .mat/.tif files;
%  vn  : name of variable in the netcdf file;
%  onm : a user-assigned name for the outputted CHIRP/CHIRPS files as character;
%  xl  : longitude (in the range of [-180 180]) of the west boundary (xl can
%        have an optional second element to represent the west boundary coordinate
%        in the unit of the output coordinate system);
%  xr  : similar to xl but for the east boundary;
%  yb  : latitude (in the range of [-50 50]) of the south boundary (yb can have
%        an optional second element to represent the south boundary coordinate
%        in the unit of the output coordinate system);
%  yt  : similar to yb but for the north boundary;

% pflg: parallel flag (false - default, squential; true - parallel);
% ors : output coordinate system (e.g. EPSG:102009);
%  rs : x and y resolution of the outputted precipitation.

%% Output
% Ofn: name list of the output .mat or .tif files stores in opth;

% onm-yyyymmddhh.mat / onm-yyyymmddhh.tif - .mat or .tif files of the processed
%  precipitation in opth;
% Grid-yyyy.mat - the lat/lon grids of the .mat file.

%% Additional note
% 1)If any optional functionality is required, please make sure to have GDAL
%   installed;
% 2)The no-data value of CHIRP/CHIRPS are preseved in the .tif files; and
% 3)Require matV2tif.m and doy2date.m.

function Ofn=GLM_process(fname,wkpth,opth,vn,onm,xl,xr,yb,yt,varargin)
%% Check the inputs
narginchk(9,12);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'fname',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'fname'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'vn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vn'));
addRequired(ips,'onm',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'onm'));
msg=sprintf('Size of xl, xr, yb, yt must be 1 or 2');
addRequired(ips,'xl',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'xr',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yb',@(x) assert(~isempty(x) & length(x)<3,msg));
addRequired(ips,'yt',@(x) assert(~isempty(x) & length(x)<3,msg));

addOptional(ips,'ors','wgs84',@(x) validateattributes(x,{'char'},{},mfilename,'ors'));
addOptional(ips,'rs',[],@(x) validateattributes(x,{'double'},{},mfilename,'rs'));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,fname,wkpth,opth,vn,onm,xl,xr,yb,yt,varargin{:});
ors=ips.Results.ors;
rs=ips.Results.rs;
pflg=ips.Results.pflg;
clear ips msg varargin

%% Lat/lon grids and other info of CHIRP/CHIRPS
Lat=double(ncread(fname,'lat'));
Lon=double(ncread(fname,'lon'));
rs_lon=360/length(Lon);
rs_lat=180/length(Lat);
Lon=0:rs_lon:360;
Lat=90:-rs_lat:-90;

[~,ys,~]=fileparts(fname);
ys=cell2mat(regexp(ys,'_(\d{4})_GLEAM','tokens','once'));

T=double(ncread(fname,'time'));
ndv=ncinfo(fname,vn);
ndv=double(ndv.FillValue); % no-data value

%% Index of interested domain
if xl(1)<0 % Convert longitude to the range of [0 360];
  xl(1)=xl(1)+360;
end
if xr(1)<0
  xr(1)=xr(1)+360;
end

cl=find(xl(1)-Lon>=0,1,'last'); % left column
cr=find(xr(1)-Lon<=0,1,'first')-1; % right column
rt=find(yt(1)-Lat<=0,1,'last'); % top row
rb=find(yb(1)-Lat>=0,1,'first')-1; % bottom row

Grid.lat=Lat(rt:rb)'-rs_lat/2;
Grid.lon=Lon(cl:cr)+rs_lon/2;
save(fullfile(opth,sprintf('Grid-%s.mat',ys)),'Grid');

xll=(cl-1)*rs_lon; % longitude of lower left corner
xll(xll>180)=xll(xll>180)-360;
yll=90-rb*rs_lat; % latitude of lower left corner

%% Parameters for resample/project/crop the file
pr1=sprintf('-t_srs %s',ors); % Project
pr2=[];
if ~isempty(rs) % Resample
  pr2=sprintf('-tr %i %i',rs(1),rs(2));
end
if length(xl)==2 && length(xr)==2 && length(yt)==2 && length(yb)==2
  pr3=sprintf('-te %i %i %i %i',xl(2),yb(2),xr(2),yt(2));
elseif length(xl)==1 && length(xr)==1 && length(yt)==1 && length(yb)==1
  pr3=[]; % Crop projected image
else
  error('sizes of xl, yb, xr, and yt must be the same and equal to 1 or 2');
end

%% Read, crop, and resample/project/crop the record
Ofn=cell(length(T),1);
switch pflg
  case true
    parfor d=1:length(T)
      Ofn{d}=GLM_process_sub(fname,vn,[1,1,d],[length(Lat)-1 length(Lon)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,onm,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
    end

  case false
    for d=1:length(T)
      Ofn{d}=GLM_process_sub(fname,vn,[1,1,d],[length(Lat)-1 length(Lon)-1 1],rt,rb,cl,cr,...
          ndv,d,ys,opth,onm,pr1,pr2,pr3,xll,yll,rs_lat,wkpth);
      fprintf('Day %s%03i processed\n',ys,d); 
    end
end

[~,~,fex]=fileparts(Ofn{1});
if strcmp(fex,'.tif')
  delete(fullfile(opth,sprintf('Grid-%s.mat',ys))); % Delete the Grid-yyyy.mat
end
end

function ofn=GLM_process_sub(fname,vn,sid,cts,rt,rb,cl,cr,ndv,d,ys,opth,onm,pr1,pr2,pr3,...
    xll,yll,rso,wkpth)
et=ncread(fname,vn,sid,cts);
et=[et(:,721:1440) et(:,1:720)]; % convert from [-180 180] to [0 360]
et=et(rt:rb,cl:cr); % crop

ds=datestr(doy2date(d,str2double(ys)),'yyyymmdd');
nm=sprintf('%s-%s',onm,ds);
ofn=fullfile(opth,sprintf('%s.mat',nm));
save(ofn,'et');

if ~contains(pr1,'wgs84') || ~isempty(pr2) || ~isempty(pr3)
  if system('gdalinfo --version')~=0
    error('GDAL is not detected. Please install GDAL to evoke the optional functionalities\n');
  end
  prm=sprintf('%s %s %s',pr1,pr2,pr3);
  prm=strrep(prm,'  ',' ');

  tfn=fullfile(wkpth,sprintf('%s.tif',nm));
  et(isnan(et))=ndv;
  matV2tif(tfn,et,xll,yll,rso,ndv,'wgs84',wkpth);
  delete(ofn); % delete the .mat file

  fun='gdalwarp -overwrite -of GTiff -r bilinear -q'; % GDAL function
  ofn=fullfile(opth,sprintf('%s.tif',nm));
  system(sprintf('%s %s "%s" "%s"',fun,prm,tfn,ofn)); % resample/project/crop
  delete(tfn);
end
end

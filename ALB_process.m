% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/9/2019

%% Functionality
% This function includes several functionalities:
%  1)read all MCD43A3 albedo tiles for the study area (resizeimg.m and matV2tif.m);
%  2)mosaic all albedo tiles together (MODISimg.m);
%  3)project, crop, resample, and output the mosaiced image (MODISimg.m and doy2date.m);

%% Input
%  afl : full name list of the MCD43A3 files;
%  arf : field name of variable (e.g. 'albedo_BSA_nir');
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: path to store the output images (e.g. C:\...\Albedo\);
% Bound: boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;
%  ors : output coordinate system (e.g. 'EPSG:102012');

% thr : a percentage of nodata value for a block where block with this percent
%       of nodata value is assigned the nodata value;
% ndv : no-data value assigned to the output images.
% pflg: parallel flag (true - default, use parallel channel; false - squential);

%% Output
% alb: processed albedo image read as Matlab variable of the study area;

% Output is the processed albedo image stored in oupth. The naming convention
%  is [arf(8:end) yyyymmdd.tif] (arf(8:end) can be, e.g., BSA_vis, WSA_nir).

%% Additional Note
% Require MODISimg.m, doy2date.m, resizeimg.m, and matV2tif.m.

function alb=ALB_process(afl,arf,wkpth,oupth,Bound,ors,varargin)
%% Check the input
narginchk(6,9);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 6 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'afl',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'afl',1));
addRequired(ips,'arf',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'arf',2));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth',3));
addRequired(ips,'oupth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'oupth',4));
addRequired(ips,'Bound',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Bound',5));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors',6));

addOptional(ips,'thr',.5,@(x) validateattributes(x,{'double'},{'scalar','>=',0,'<=',1},mfilename,'thr',7));
addOptional(ips,'ndv',-999,@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv',8));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg',9));

parse(ips,afl,arf,wkpth,oupth,Bound,ors,varargin{:});
thr=ips.Results.thr;
ndv=ips.Results.ndv;
pflg=ips.Results.pflg;
clear ips

%% Properties of input records
hif=hdfinfo(afl(1,:));
hif=hif.Vgroup.Vgroup(1).SDS;
vrl=arrayfun(@(x) x.Name,hif,'UniformOutput',false);
id=strcmp(vrl,arf);
ndv_o=double(hif(id).Attributes(4).Value); % no-data value
scf=double(hif(id).Attributes(5).Value); % multiplicative factor

hif=hdfinfo(afl(1,:),'EOS');
hif=hif.Grid;
xli=hif.UpperLeft(1);
xri=hif.LowerRight(1);
yti=hif.UpperLeft(2);
ybi=hif.LowerRight(2);
rxi=(xri-xli)/hif.Columns;
ryi=(yti-ybi)/hif.Rows;
kx=max(1,fix(Bound(3,1)/rxi));
ky=max(1,fix(Bound(3,2)/ryi));

%% Upscale the image
Tfl=[];
idn=fullfile(wkpth,['idA_r' num2str(Bound(3,1),'%i') '.mat']);
ors_i='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" ';

tfn=ALB_process_sub(afl(1,:),arf,ndv_o,scf,kx,ky,ndv,idn,thr,wkpth,1,rxi,ors_i);
Tfl=[Tfl;tfn];
switch pflg
  case true
    parfor n=2:size(afl,1)
      tfn=ALB_process_sub(afl(n,:),arf,ndv_o,scf,kx,ky,ndv,idn,thr,wkpth,n,rxi,ors_i);
      Tfl=[Tfl;tfn];
    end

  case false
    for n=2:size(afl,1)
      tfn=ALB_process_sub(afl(n,:),arf,ndv_o,scf,kx,ky,ndv,idn,thr,wkpth,n,rxi,ors_i);
      Tfl=[Tfl;tfn];
    end
end

%% Process the upsacled image
[~,ds,~]=fileparts(afl(1,:));
[ds,~]=regexp(ds,'.A(\d{4})(\d{3}).','tokens','once');
ys=ds{1};
ds=doy2date(str2double(ds{2}),str2double(ys));
ds=datestr(ds,'yyyymmdd');

% Image processing (read, mosaic, project, crop, resample)
oun=fullfile(oupth,[arf(8:end) ds '.tif']);
alb=MODISimg(Tfl,wkpth,oun,Bound,ors,'','','bilinear',pflg);
end

function tfn=ALB_process_sub(afn,arf,ndv_o,scf,kx,ky,ndv,idn,thr,wkpth,n,rxi,ors_i)
% Upscale image
alb=double(hdfread(afn,arf));
alb(alb==ndv_o)=NaN;
alb=alb*scf;
if max(kx,ky)>1
  alb=resizeimg(alb,ndv,kx,ky,idn,thr);
end

% Write the upscaled image to the working directory
hif=hdfinfo(afn,'EOS');
hif=hif.Grid;
xli=hif.UpperLeft(1);
ybi=hif.LowerRight(2);
[~,nm,~]=fileparts(afn);
nm=cell2mat(regexp(nm,'A(?<year>\d+)(?<day>\d{3})','match'));
tfn=sprintf('%s%02i.tif',[wkpth nm '_p'],n);

matV2tif(tfn,alb,xli,ybi,kx*rxi,ndv,ors_i,wkpth);
end

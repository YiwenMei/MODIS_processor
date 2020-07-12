% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 4/21/2020

%% Functionality
% This function is used to process the SNODAS product (NOHRSC 2004). Its main
%  functionalities are
%   1)untar the eight variables for one time step stored in the daily file;
%   2)unzip selected variable(s) from the eight; and
%   3)convert the selected variable(s) in .dat to .tif format.
%  The extracted variable(s) are outputted as .tif files to the output directory.
%  Multiple optional functionalities listed below may be included
%   4)crop the variable(s);
%   5)resample the variable(s); and/or
%   6)project the variable(s) to another coordinate system.
%  The processed files are outputted as .tif in the output directory.

%% Input
%  tfn : full name of the output geotiff file;
% wkpth: working directory of the code;
% opth : output directory of the processed files;
%  vid : id of the record of interest;
%  vno : output names of the record of interest;
%  GIf : boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;

% ors: coordinate system of the files;
% itm: interpolation method accepted by gdalwarp (e.g. 'bilinear');

%% Output
% The output files are located in opth with the name "SNODAS-vno-yyyymmdd.tif".

function SDAS_process(tfn,wkpth,opth,vid,vno,GIf,varargin)
%% Check the input
narginchk(6,8);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'tfn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'tfn'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
if length(vid)~=length(vno)
  error('Number of element in vid must be correspond to those in vno');
end
addRequired(ips,'vid',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'vid'));
addRequired(ips,'vno',@(x) validateattributes(x,{'cell'},{'vector'},mfilename,'vno'));
addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3 2]},mfilename,'GIf'));

pr1='-a_srs "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"';
addOptional(ips,'ors',pr1,@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));
addOptional(ips,'itm','bilinear',@(x) any(strcmp(x,{'near','bilinear','cubic','cubicspline',...
    'lanczos','average','mode','max','min','med','q1','q3'})));

parse(ips,tfn,wkpth,opth,vid,vno,GIf,varargin{:});
ors=ips.Results.ors;
itm=ips.Results.itm;
clear ips varargin

%% File header
L1='ENVI';
L2={'samples=6935','samples=8192'};
L3={'lines=3351','lines=4096'};
L4='bands=1';
L5='header offset=0';
L6='file type=ENVI Standard';
L7='data type=2';
L8='interleave=bsq';
L9='byte order=1';

[~,ds,~]=fileparts(tfn);
ty=cell2mat(regexp(ds,'_(\w+)_','tokens','once'));
switch ty
  case 'masked'
    ty=1;
  case 'unmasked'
    ty=2;
  otherwise
    error('Please input the daily ".tar" file');
end
ds=cell2mat(regexp(ds,'_(\d{8})','tokens','once'));
if datenum(ds,'yyyymmdd')<=datenum('20131001','yyyymmdd')
  yf=2;
else
  yf=1;
end

%% GDAL functions
fun1='gdal_translate -q -of GTiff';
pr2='-a_nodata -9999';
pr3={'-a_ullr -124.73333333 52.87500000 -66.94166667 24.95000000',...
    '-a_ullr -130.516666666661 58.2333333333310 -62.2499999999975 24.0999999999990';...
    '-a_ullr -124.73375000 52.87458333 -66.94208333 24.87458333',...
    '-a_ullr -130.517083333328 58.2329166666644 -62.2504166666637 24.0995833333324'};
prm1=sprintf('%s %s %s',pr1,pr2,pr3{yf,ty});

fun2='gdalwarp -of GTiff -overwrite -q';
pr1=sprintf('-t_srs %s',ors);
pr3=sprintf('-r %s',itm);
if ~isnan(GIf(3,1)) && ~isnan(GIf(1,1))
  pr2=sprintf('-te %i %i %i %i',GIf(1,1),GIf(2,2),GIf(2,1),GIf(1,2));
  pr4=sprintf('-tr %i %i',GIf(3,1),GIf(3,2));
  prm2=sprintf('%s %s %s %s',pr1,pr2,pr3,pr4);
elseif ~isnan(GIf(3,1)) && isnan(GIf(1,1))
  pr4=sprintf('-tr %i %i',GIf(3,1),GIf(3,2));
  prm2=sprintf('%s %s %s',pr1,pr3,pr4);
elseif isnan(GIf(3,1)) && ~isnan(GIf(1,1))
  pr2=sprintf('-te %i %i %i %i',GIf(1,1),GIf(2,2),GIf(2,1),GIf(1,2));
  prm2=sprintf('%s %s %s',pr1,pr2,pr3);
else
  prm2=sprintf('%s %s %s',pr1,pr3);
end

%% Untar the daily file and unzip the variable file
tfl=untar(tfn,wkpth);
k=cellfun(@(X) contains(X,'.dat.gz'),tfl);
cellfun(@delete,tfl(~k)); % delete the .txt.gz
tfl=tfl(k)';
for v=1:length(vid)
  fn=cell2mat(gunzip(tfl{vid(v)}));
  [~,ds,~]=fileparts(fn);
  ds=cell2mat(regexp(ds,'TNATS(\d{10})[HD]','tokens','once'));

%% Create header file
  afn=strrep(fn,'.dat','.hdr');
  fid=fopen(afn,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',L1,L2{ty},L3{ty},L4,L5,L6,L7,L8,L9);
  dlmwrite(afn,'','-append');
  fclose(fid);

%% Convert the .dat to .tif
  ifn=fullfile(wkpth,sprintf('SNODAS-%s-%s.tif',vno{v},ds));
  cmstr=sprintf('%s %s %s %s',fun1,prm1,fn,ifn);
  system(cmstr);
  delete(fn); % delete the .dat and .hdr
  delete(afn);

%% Crop, project, and resample
  ofn=fullfile(opth,sprintf('SNODAS-%s-%s.tif',vno{v},ds));
  if strcmp(ors,pr1) && all(isnan(GIf))
    movefile(ifn,ofn);
  else
    cmstr=sprintf('%s %s %s %s',fun2,prm2,ifn,ofn);
    system(cmstr);
    delete(ifn); % delete the converted .tif
  end

  fprintf('Variable "%s" of step "%s" processed\n',vno{v},ds);
end
cellfun(@delete,tfl); % delete the .dat.gz
end

% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 8/5/2019

%% Functionality
% This function includes several functionalities:
%  1)read all images for the study area;
%  2)mosaic all images together;
%  3)reproject, crop and resample the mosaiced image;

%% Input
%  imL : full name list of the MODIS hdf files (e.g.
%        C:\...\MCD12Q1.A2007001.h23v05.006.2018054173442.hdf;
%        C:\...\MOD11A1.A2007001.h23v05.006.2015307155340.hdf);
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
%  oun : full name of output geotiff image (e.g. C:\...\XXX.tif; set it to "[]"
%        if no need to output the image);
%  GIf : boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;
%  ors : coordinate system of output image (e.g. 'EPSG:102012');

%  rn : record name of the MODIS image (e.g. 'MODIS_Grid_Daily_1km_LST' for MOD11A1;
%       'MCD12Q1' for MCD12Q1);
% vrf : field of interest of the images (e.g. 'LC_Type3' for MCD12Q1; 'Emis_31'
%       for MYD11A1);
% itm : interpolation method accepted by gdalwarp (e.g. 'bilinear');
% pflg: parallel flag (true - default, use parallel channel; false - squential).

%% Output
% imout: output image (the image file is also outputted as "oun").

%% Additional note
% 1)Need to have GDAL installed to run the code;
% 2)No-data value of the original image preserved in "oun".

function imout=MODISimg(imL,wkpth,GIf,ors,varargin)
%% Check the input
narginchk(4,9);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'iml',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'iml'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3 2]},mfilename,'GIf'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));

addOptional(ips,'oun','',@(x) validateattributes(x,{'char'},{},mfilename,'oun'));
addOptional(ips,'rn','',@(x) validateattributes(x,{'char'},{},mfilename,'rn'));
addOptional(ips,'vrf','',@(x) validateattributes(x,{'char'},{},mfilename,'vrf'));
addOptional(ips,'itm','bilinear',@(x) any(strcmp(x,{'near','bilinear','cubic','cubicspline',...
    'lanczos','average','mode','max','min','med','q1','q3'})));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg'));

parse(ips,imL,wkpth,GIf,ors,varargin{:});
oun=ips.Results.oun;
rn=ips.Results.rn;
vrf=ips.Results.vrf;
itm=ips.Results.itm;
pflg=ips.Results.pflg;
clear ips varargin

%% Convert original to geotiff
[~,nm,fex]=fileparts(imL(1,:));
if ~strcmp(fex,'.tif')
  if isempty(rn) || isempty(vrf)
    error('rn and vrf cannot be empty');
  else
    ds=cell2mat(regexp(nm,'.A(\d{7}).h','tokens','once'));

    fun='gdal_translate -of GTiff -q '; % GDAL function
    pr1=sprintf('-r %s',itm);
    switch pflg
      case true
        parfor n=1:size(imL,1)
          [a,~]=regexp(imL(n,:),'h+\d{2}v+\d{2}','match');
          a=cell2mat(a);
          inv=sprintf('HDF4_EOS:EOS_GRID:"%s":%s:%s',imL(n,:),rn,vrf); % Full name with record
          im1=fullfile(wkpth,sprintf('%s-%s-%s.tif',vrf,ds,a)); % name and field name

          system(sprintf('%s %s %s "%s"',fun,pr1,inv,im1));
        end

      case false
        for n=1:size(imL,1)
          [a,~]=regexp(imL(n,:),'h+\d{2}v+\d{2}','match');
          a=cell2mat(a);
          inv=sprintf('HDF4_EOS:EOS_GRID:"%s":%s:%s',imL(n,:),rn,vrf); % Full name with record
          im1=fullfile(wkpth,sprintf('%s-%s-%s.tif',vrf,ds,a)); % name and field name

          system(sprintf('%s %s %s "%s"',fun,pr1,inv,im1));
        end
    end
    inv=fullfile(wkpth,sprintf('%s-%s-h*v*.tif',vrf,ds));
  end
  
else
  inv=imL;
  wkpth=fileparts(inv(1,:));
  nm=strsplit(nm,'-');
  vrf=nm{1};
  ds=nm{2};
  inv=fullfile(wkpth,sprintf('%s-%s-h*v*.tif',vrf,ds));
end

%% Build virtual mosaiced image
fun='gdalbuildvrt -overwrite -q';
im1=fullfile(wkpth,sprintf('%s-%s.vrt',vrf,ds));

cmstr=sprintf('%s %s %s',fun,im1,inv);
system(cmstr);

%% Project image
fun='gdalwarp -of GTiff -overwrite -q';
pr1=sprintf('-t_srs %s',ors);
pr2=sprintf('-te %i %i %i %i',GIf(1,1),GIf(2,2),GIf(2,1),GIf(1,2));
pr3=sprintf('-r %s',itm);
if ~isnan(GIf(3,1)) && ~isnan(GIf(3,2))
  pr4=sprintf('-tr %i %i',GIf(3,1),GIf(3,2));
  prm=sprintf('%s %s %s %s',pr1,pr2,pr3,pr4);
else
  prm=sprintf('%s %s %s',pr1,pr2,pr3);
end
ouv=fullfile(wkpth,sprintf('%s-%s.tif',vrf,ds));

cmstr=sprintf('%s %s "%s" "%s"',fun,prm,im1,ouv);
system(cmstr);
delete(inv);
delete(im1);

%% Move the processed image to ouv
imout=double(imread(ouv));
if isempty(oun)
  delete(ouv);
else
  movefile(ouv,oun);
end
end

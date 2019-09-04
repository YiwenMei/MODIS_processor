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
% Bound: boundary latitude and longitude and resolution of output images (it
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

function imout=MODISimg(imL,wkpth,oun,Bound,ors,varargin)
%% Check the input
narginchk(5,9);
ips=inputParser;
ips.FunctionName=mfilename;
fprintf('%s received 7 required and %d optional inputs\n',mfilename,length(varargin));

addRequired(ips,'iml',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'iml',1));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth',2));
addRequired(ips,'oun',@(x) validateattributes(x,{'char'},{},mfilename,'wkpth',3));
addRequired(ips,'Bound',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'Bound',4));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors',5));

addOptional(ips,'rn','',@(x) validateattributes(x,{'char'},{},mfilename,'rn',6));
addOptional(ips,'vrf','',@(x) validateattributes(x,{'char'},{},mfilename,'vrf',7));
addOptional(ips,'itm','bilinear',@(x) any(strcmp(x,{'near','bilinear','cubic','cubicspline',...
    'lanczos','average','mode','max','min','med','q1','q3'})));
addOptional(ips,'pflg',false,@(x) validateattributes(x,{'logical'},{'nonempty'},mfilename,'pflg',9));

parse(ips,imL,wkpth,oun,Bound,ors,varargin{:});
rn=ips.Results.rn;
vrf=ips.Results.vrf;
itm=ips.Results.itm;
pflg=ips.Results.pflg;
clear ips

%% Convert original to geotiff
[~,nm,fex]=fileparts(imL(1,:));
if ~strcmp(fex,'.tif')
  if isempty(rn) || isempty(vrf)
    error('rn and vrf cannot be empty');
  else
    nm=regexp(nm,'A(?<year>\d{4})(?<day>\d{3})','match');
    nm=cell2mat(nm);

    fun='gdal_translate -of GTiff '; % GDAL function
    pr1=['-r ' itm ' '];
    switch pflg
      case true
        parfor n=1:size(imL,1)
          inv=['HDF4_EOS:EOS_GRID:"' imL(n,:) '":' rn ':' vrf]; % Full name with record
          im1=fullfile(wkpth,[nm '_p' num2str(n,'%02i') '.tif']); % name and field name

          system([fun pr1 inv ' "' im1 '"']);
        end

      case false
        for n=1:size(imL,1)
          inv=['HDF4_EOS:EOS_GRID:"' imL(n,:) '":' rn ':' vrf]; % Full name with record
          im1=fullfile(wkpth,[nm '_p' num2str(n,'%02i') '.tif']); % name and field name

          system([fun pr1 inv ' "' im1 '"']);
        end
    end
    inv=fullfile(wkpth,[nm '_p*.tif']);
  end
  
else
  inv=imL;
  fpth=fileparts(inv(1,:));
  nm=cell2mat(regexp(nm,'A(?<year>\d{4})(?<day>\d{3})','match'));
  inv=fullfile(fpth,[nm '_p*.tif']);
end

%% Build virtual mosaiced image
fun='gdalbuildvrt -overwrite ';
im1=fullfile(wkpth,[nm '.vrt']);

system([fun '"' im1 '" ' inv]); % On linux
% system([fun '"' im1 '" "' inv '"']); % On windows

%% Project image
fun='gdalwarp -of GTiff -overwrite ';
pr1=['-t_srs ' ors ' '];
pr2=sprintf('-te %i %i %i %i ',Bound(1,1),Bound(2,2),Bound(2,1),Bound(1,2));
pr3=['-r ' itm ' '];
if ~isnan(Bound(3,1)) && ~isnan(Bound(3,2))
  pr4=sprintf('-tr %i %i ',Bound(3,1),Bound(3,2));
else
  pr4=[];
end
prm=[pr1 pr2 pr3 pr4];
ouv=fullfile(wkpth,[nm '.tif']);

system([fun prm '"' im1 '" "' ouv '"']);
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

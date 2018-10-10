% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/05/2018

%% Functionality
% This function includes several functionalities:
%  1)read all images for the study area;
%  2)mosaic all images together;
%  3)reproject, crop and resample the mosaiced image;

%% Input
%  imL : full name list of the MODIS hdf files (e.g.
%        C:\...\MCD12Q1.A2007001.h23v05.006.2018054173442.hdf;
%        C:\...\MOD11A1.A2007001.h23v05.006.2015307155340.hdf);
%  rn  : record name of the MODIS image (e.g. 'MODIS_Grid_Daily_1km_LST' for
%        MOD11A1; 'MCD12Q1' for MCD12Q1);
%  vrf : field of interest of the images (e.g. 'LC_Type3' for MCD12Q1; 'Emis_31'
%        for MYD11A1);
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
%  oun : full name of output geotiff image (e.g. C:\...\XXX.tif; set it to "[]"
%        if no need to output the image);
% rx/ry: resolution of output images;
% xl/xr: left/right side coordinate;
% yb/yt: same as xl/xr but for the bottom/top side's;
%  ors : coordinate system of output image (e.g. 'EPSG:102012');
%  itm : interpolation method accepted in GDAL (e.g. 'bilinear').

%% Output
% imout: output image (the image file is also outputted as "oun").

%% Additional note
%   1)Need to have GDAL installed so as to run the code;
%   2)No-data value of the original image preserved in "oun".

function imout=MODISimg(imL,rn,vrf,wkpth,oun,xl,xr,rx,yb,yt,ry,ors,itm)
%% Convert original to geotiff
[~,~,fex]=fileparts(imL(1,:));

if ~strcmp(fex,'.tif')
  fun='gdal_translate -of GTiff '; % GDAL function
  pr1=['-r ' itm ' '];
  for n=1:size(imL,1)
    inv=['HDF4_EOS:EOS_GRID:"' imL(n,:) '":' rn ':' vrf]; % Full name with record
    im1=fullfile(wkpth,['p' num2str(n,'%02i') '.tif']);   % name and field name

    system([fun pr1 inv ' "' im1 '"']);
  end

  inv=fullfile(wkpth,'p*.tif');
else
  inv=imL;
end
%% Build virtual mosaiced image
fun='gdalbuildvrt -overwrite ';
im1=fullfile(wkpth,'p.vrt');

system([fun '"' im1 '" ' inv]); % On linux
% system([fun '"' im1 '" "' inv '"']); % On windows

%% Project image
fun='gdalwarp -of GTiff -overwrite -r bilinear ';
pr1=['-t_srs ' ors ' '];
pr2=sprintf('-te %i %i %i %i ',xl,yb,xr,yt);
if ~isempty(rx) && ~isempty(ry)
  pr3=sprintf('-tr %i %i ',rx,ry);
else
  pr3=[];
end
par=[pr1 pr2 pr3];
ouv=fullfile(wkpth,'imout.tif');

system([fun par '"' im1 '" "' ouv '"']);
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

% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/07/2018

%% Functionality
% This function contains several functionalities:
%  1)read all MCD12Q1 LC Type 3 images for the study area (MODISimg.m);
%  2)mosaic all LC images together (MODISimg.m);
%  3)reproject, crop and resample the mosaiced LC image (MODISimg.m);
%  4)calculate the fraction of a land cover class for a block (LC_Frac.m);
%  5)output the land cover and land cover fraction images.

%% Input
% lcfl: full name list of the land cover files;
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: output directory for the images (e.g. C:\...\LCFrac\);
% rx/ry: resolution of output images;
% xl/xr: left/right side coordinate;
% yb/yt: same as xl/xr but for the bottom/top side's;
%  ors : output coordinate system (e.g. 'EPSG:102012');
%  thr : a percentage of nodata value for a block where block with this percent
%        of nodata value is assigned the nodata value;
%  ndv : no-data value assigned to the output images.

%% Output
% imout: mosaiced land cover image of the study area;
% Land cover fraction images are stored in oupth with the naming convention as
%  LCXX_yyyy.tif where XX represents the ID of a land cover class.

%% Additional Note
%  This function processes the MCD12Q1 land cover type 3 record. ID of the land
%   cover classes are:
%   00 - water body,          01 - Grassland,            02 - Shrublands,
%   03 - Broadleaf Croplands, 04 - Savannas,             05 - Evergreen Broadleaf,
%   06 - Deciduous Broadleaf, 07 - Evergreen Needleleaf, 08 - Deciduous Needleleaf,
%   09 - Unvegetated,         10 - Urban and Built-up Lands.

function imout=LC_process(lcfl,wkpth,oupth,xl,xr,rx,yb,yt,ry,ors,thr,ndv)
%% Properties of input records
hif=hdfinfo(lcfl(1,:));
rn=hif.Vgroup.Name; % Name of the record
hif=hif.Vgroup.Vgroup(1).SDS(3); % LC type 3
lcc=hif.Attributes(2:12); % Land cover classes
vrf=hif.Name; % Name of the field
ndv_o=double(hif.Attributes(13).Value); % no-data-value of LC

%% Image processing (read, mosaic, project, crop, resample)
imout=MODISimg(lcfl,rn,vrf,wkpth,[],xl,xr,[],yb,yt,[],ors,'near');
imout(imout==ndv_o)=NaN;

%% Calculate land cover fraction
[~,ys,~]=fileparts(lcfl(1,:));
ys=ys(10:13);
imo=[wkpth 'p.asc'];

rsx_in=(xr-xl)/size(imout,2);
rsy_in=(yt-yb)/size(imout,1);
for v=1:length(lcc)
  LC=LC_Frac(imout,[xl xr rsx_in],[xl xr rx],[yt yb rsy_in],[yt yb ry],...
      fullfile(wkpth,sprintf('id%irs.mat',rx)),double(lcc(v).Value),thr,ndv);

% Output the land cover fraction images
  fid=fopen(imo,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(LC,2))],...
      ['nrows ' num2str(size(LC,1))],['xllcorner ' num2str(xl,12)],...
      ['yllcorner ' num2str(yb,12)],['cellsize ' num2str(rx,12)],...
      ['NODATA_value ' num2str(ndv)]);
  dlmwrite(imo,LC,'delimiter',' ','-append');
  fclose(fid);

  fun='gdal_translate -of GTiff -r bilinear '; % GDAL function
  pr1=['-a_srs ' ors ' '];
  IMo=fullfile(oupth,sprintf('LC%02i_%s.tif',lcc(v).Value,ys));

  system([fun pr1 imo ' "' IMo '"']);
  delete(imo);
end
end

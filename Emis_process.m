% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/05/2018

%% Functionality
% This function includes several functionalities:
%  1)read all the MOD11A1 band 31 & 32 emissivity images for the study area
%    (MODISimg.m);
%  2)mosaic all the emissivity images together (MODISimg.m);
%  3)project, crop and resample the mosaiced emissivity image (MODISimg.m);
%  4)calculate the mean of band 31 and band 32 emissivity image;
%  5)calculate the mean of MOD and MYD -based band-wise mean emissivity image
%    attained from step 4.

%% Input
% Emfl : a 2-by-1 cell array with cell 1/2 stores a list of full names of the
%        MOD11A1/MYD11A1 records for the location of interest for a day;
% wkpth: working directory for the code (e.g. C:\...\wkdir\);
% oupth: output directory for the images (e.g. C:\...\Emis\2001\);
% rx/ry: resolution of output emissivity with same unit required by "coor";
% xl/xr: left/right corner coordinate with same unit required by "coor";
% yb/yt: same as xl/xr but for the bottom and top corner coordinate;
% coor : coordinate system of output images;
%  ndv : no-data value assigned to the output images.

%% Output
% Output images are stored in oupth as EMSyyyymmdd.tif.

function Emis_process(Emfl,wkpth,oupth,xl,xr,rx,yb,yt,ry,coor,ndv)
emfl=cell2mat(Emfl); % Check whether the list is empty
if ~isempty(emfl)
%% Properties of input records
  hif=hdfinfo(emfl(1,:));
  rn=hif.Vgroup.Name; % Name of the record
  vrf=[hif.Vgroup.Vgroup(1).SDS(9).Name;hif.Vgroup.Vgroup(1).SDS(10).Name];
  hif=hif.Vgroup.Vgroup(1).SDS(9);
  scf=double(hif.Attributes(6).Value); % scale factor
  ofs=double(hif.Attributes(7).Value); % additive offset
  ndv_o=double(hif.Attributes(4).Value); % no-data value of VI

%% Process the emissivity data
  IMO=nan(round(yt-yb)/ry,round(xr-xl)/rx);
  for t=1:length(Emfl)
    IMo=nan(round(yt-yb)/ry,round(xr-xl)/rx);

    for b=1:size(vrf,1)
      if ~isempty(Emfl{t})
% Image processing (read, mosaic, project, crop, resample)
        imo=MODISimg(Emfl{t},rn,vrf(b,:),wkpth,[],xl,xr,rx,yb,yt,ry,coor,'bilinear');
        imo(imo==ndv_o)=NaN;
        imo=imo*scf+ofs;
      else
        imo=nan(round(yt-yb)/ry,round(xr-xl)/rx);
      end

% Mean of emissivity
      IMo(:,:,b)=imo;
    end
    IMo=nanmean(IMo,3); % Mean of two bands
    IMO(:,:,t)=IMo;
  end
  IMO=nanmean(IMO,3); % Mean of MOD and MYD
  IMO(isnan(IMO))=ndv;

%% Output the processed emissivity
  imo=[wkpth 'p.asc'];
  [~,ns,~]=fileparts(emfl(1,:));
  ys=ns(10:13);
  ds=ns(14:16);
  ds=doy2date(str2double(ds),str2double(ys));
  ds=datestr(ds,'yyyymmdd');

  fid=fopen(imo,'w');
  fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(IMO,2))],['nrows '...
      num2str(size(IMO,1))],['xllcorner ' num2str(xl,12)],['yllcorner '...
      num2str(yb,12)],['cellsize ' num2str(rx,12)],['NODATA_value ' num2str(ndv)]);
  dlmwrite(imo,IMO,'delimiter',' ','-append');
  fclose(fid);

  fun='gdal_translate -of GTiff -r bilinear '; % GDAL function
  pr1=['-a_srs ' coor ' '];
  IMo=fullfile(oupth,['EMS' ds '.tif']);

  system([fun pr1 '"' imo '" "' IMo '"']);
  delete(imo);
end
end

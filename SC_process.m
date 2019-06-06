% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/05/2018

%% Functionality
% This function includes several functionalities:
%  1)read and average the MOD10A1 and MYD10A1 NDSI Snow Cover images;
%  2)output the mean SC to the temp folder;
%  3)read and mosaic all the SC images together (MODISimg.m);
%  4)project, crop and resample the mosaiced SC image (MODISimg.m);

%% Input
% LTfl : a 2-by-1 cell array with cell 1/2 stores a list of full names of the
%        MOD11A1/MYD11A1 records for the location of interest for a day;
% wkpth: working directory, the temp folder, for the code (e.g. C:\...\wkdir\);
% opth : output directory for the images (e.g. C:\...\Emis\2001\);
% rx/ry: resolution of output images;
% xl/xr: left/right side coordinate;
% yb/yt: same as xl/xr but for the bottom/top side's;
%  ors : output coordinate system (e.g. 'EPSG:102012');
%  ndv : no-data value assigned to the output images.

%% Output
% imo: final SC data;
% Output images are stored in oupth as SCyyyymmdd.tif.

function imo=SC_process(LTfl,wkpth,opth,xl,xr,rx,yb,yt,ry,ors,ndv)
%% Check the input
switch nargin
  case {1:10}; error('Not enough number of arguments');
  case 11
  otherwise; error('Too many number of arguments');
end

%% Properties of input records
LTfl=cell2mat(LTfl); % Check whether the list is empty
ors_o='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"';
if ~isempty(LTfl)
  hif=hdfinfo(LTfl(1,:));
  vrf=hif.Vgroup.Vgroup(1).SDS(1).Name;
  hif=hif.Vgroup.Vgroup(1).SDS(1);
  ndv_o=double(hif.Attributes(5).Value); % no-data value of VI

%% Find unique tiles
  A=cell(size(LTfl,1),1);
  for t=1:size(LTfl,1)
    [a,~]=regexp(LTfl(t,:),'h+\d{2}v+\d{2}','match');
    A{t}=cell2mat(a);
  end
  A=unique(A);

%% Processing
  Tfl=[];
  LTfl=mat2cell(LTfl,ones(size(LTfl,1),1),size(LTfl,2));
  parfor n=1:length(A)
    k=strfind(LTfl,A{n});
    k=cellfun(@isempty,k);
    hfl=LTfl(~k);

% Average the MOD and MYD record
    LT=[];
    for t=1:length(hfl)
      lt=double(hdfread(hfl{t},vrf));
      lt(lt>=ndv_o)=NaN;
      LT=cat(3,LT,lt);
    end
    LT=nanmean(LT,3);
    LT(isnan(LT))=ndv;

% Output the averaged LST of a tile
    hif=hdfinfo(hfl{1},'eos');
    xll=hif.Grid.UpperLeft(1,1);
    yll=hif.Grid.LowerRight(1,2);
    rs_o=mean(abs(hif.Grid.UpperLeft-hif.Grid.LowerRight)./[hif.Grid.Columns hif.Grid.Rows]);
    [~,nm,~]=fileparts(hfl{1});
    nm=regexp(nm,'A(?<year>\d+)(?<day>\d{3})','match');
    nm=cell2mat(nm);
    tfn=sprintf('%s%02i.tif',[wkpth nm '_p'],n);
    matV2tif(tfn,LT,xll,yll,rs_o,ndv,ors_o,wkpth);

% Image processing (read, mosaic, project, crop, resample)
    Tfl=[Tfl;tfn];
  end
  tfl=mat2cell(Tfl,ones(size(Tfl,1),1),size(Tfl,2));
  tfl=regexp(tfl,'(?<year>\d{4})(?<day>\d{3})','match');
  ds=cellfun(@cell2mat,tfl,'UniformOutput',false);
  ds=cell2mat(unique(ds));
  ds=datestr(doy2date(str2double(ds(5:end)),str2double(ds(1:4))),'yyyymmdd');
  oun=[opth 'SC' ds '.tif'];
  imo=MODISimg(Tfl,[],[],wkpth,oun,xl,xr,rx,yb,yt,ry,ors,'bilinear');
  imo(imo==ndv)=NaN;
end
end

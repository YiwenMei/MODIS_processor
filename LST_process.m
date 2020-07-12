% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 4/28/2020

%% Functionality
% This function includes several functionalities:
%  1)read all the tiles for MOD11A1 LST_Day & LST_Night records over the study
%     area and convert them to geotiff;
%  2)mosaic all the LST geotiff files together (MODISimg.m);
%  3)project, crop and resample the mosaiced LST (MODISimg.m);
%  4)calculate the mean of MOD and MYD LST_D and LST_N (4 records) to represent
%     LST of the day.

%% Input
%  Tfl : a 2-by-1 cell array with cell 1/2 stores a list of full names of the
%        MOD11A1/MYD11A1 .hdf files for the location of interest for a time step;
% wkpth: working directory for the code (e.g., /.../wkdir/);
% opth : output directory for the images (e.g. /.../LST/);
%  GIf : boundary latitude and longitude and resolution of output images (it
%        follows [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution,
%        l/r/b/t stands for left/right/bottom/top;
%  ors : output coordinate system (e.g., 'EPSG:102012');
%  ndv : no-data value assigned to the output files.

%% Output
% Output images are stored in opth as MODIS-LST-yyyymmdd.tif.

function LST_process(Tfl,wkpth,opth,GIf,ors,ndv)
%% Check the input
narginchk(6,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Tfl',@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'Tfl'));
addRequired(ips,'wkpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wkpth'));
addRequired(ips,'opth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'opth'));
addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'GIf'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));
addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));

parse(ips,Tfl,wkpth,opth,GIf,ors,ndv);
clear ips varargin

ors_o='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"';
cid=[1 1 2 2;1 5 1 5];
Vb=[];
N=[];
for r=1:size(cid,2)
  tfl=Tfl{cid(1,r)};
  if ~isempty(tfl)
%% Convert .hdf to geotiff
    tmfl=[];
    for t=1:size(tfl,1)
      hif=hdfinfo(tfl(t,:),'eos');
      xll=hif.Grid.UpperLeft(1);
      yll=hif.Grid.LowerRight(2);
      rs=mean((hif.Grid.UpperLeft-hif.Grid.LowerRight)./[-hif.Grid.Columns hif.Grid.Rows]);

      hif=hdfinfo(tfl(t,:));
      hif=hif.Vgroup.Vgroup(1).SDS(cid(2,r));
      vn=hif.Name; % LST Day/Night
      scf=double(hif.Attributes(7).Value); % scale factor
      nv_o=double(hif.Attributes(5).Value); % no-data value of VI

      vb=double(hdfread(tfl(t,:),vn));
      vb(vb==nv_o)=NaN;
      vb=vb*scf;

      [~,nm,~]=fileparts(tfl(t,:));
      ds=cell2mat(regexp(nm,'.A(\d{7}).','tokens','once'));
      nm=cell2mat(regexp(nm,'.(h\d{2}v\d{2}).','tokens','once'));
      ofn=fullfile(wkpth,sprintf('LST-%s-%s.tif',ds,nm));
      tmfl=[tmfl;ofn];
      matV2tif(ofn,vb,xll,yll,rs,ndv,ors_o,wkpth);
    end

%% Mosaic, project, crop, and resample
    vb=MODISimg(tmfl,wkpth,GIf,ors);
    vb(vb==ndv)=NaN;
    Vb=nansum(cat(3,Vb,vb),3);
    n=~isnan(vb);
    N=sum(cat(3,N,n),3);
  end
end

%% Mean LST
vb=Vb./N;
vb(N==0)=NaN;
ds=datestr(doy2date(str2double(ds(5:7)),str2double(ds(1:4))),'yyyymmdd');
ofn=fullfile(opth,sprintf('MODIS-LST-%s.tif',ds));
matV2tif(ofn,vb,GIf(1,1),GIf(2,2),GIf(3,1),ndv,ors,wkpth);
end

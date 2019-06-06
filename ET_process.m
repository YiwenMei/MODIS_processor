function imo=ET_process(ETfl,wkpth,opth,xl,xr,rx,yb,yt,ry,ors,ndv)
%% Check the input
switch nargin
  case {1:10}; error('Not enough number of arguments');
  case 11
  otherwise; error('Too many number of arguments');
end

%% Properties of input records
etfl=cell2mat(ETfl); % Check whether the list is empty
ors_o='"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"';
if ~isempty(etfl)
  hif=hdfinfo(etfl(1,:));
  vrf=hif.Vgroup.Vgroup(1).SDS(1).Name;
  hif=hif.Vgroup.Vgroup(1).SDS(1);
  scf=double(hif.Attributes(5).Value); % scale factor
  ndv_o=double(hif.Attributes(3).Value); % no-data value of ET

%% Find unique tiles
  A=cell(size(etfl,1),1);
  for t=1:size(etfl,1)
    [a,~]=regexp(etfl(t,:),'h+\d{2}v+\d{2}','match');
    A{t}=cell2mat(a);
  end
  A=unique(A);

%% Processing
  Tfl=[];
  etfl=mat2cell(etfl,ones(size(etfl,1),1),size(etfl,2));
  parfor n=1:length(A)
    k=strfind(etfl,A{n});
    k=cellfun(@isempty,k);
    hfl=etfl(~k);

% Average the MOD and MYD record
    Et=[];
    for t=1:length(hfl)
      et=double(hdfread(hfl{t},vrf));
      et(et>max(ndv_o) | et<min(ndv_o))=NaN;
      Et=cat(3,Et,et);
    end
    Et=nanmean(Et,3)*scf;
    Et(isnan(Et))=ndv;

% Output the averaged ET of a tile
    hif=hdfinfo(hfl{1},'eos');
    xll=hif.Grid.UpperLeft(1,1);
    yll=hif.Grid.LowerRight(1,2);
    rs_o=mean(abs(hif.Grid.UpperLeft-hif.Grid.LowerRight)./[hif.Grid.Columns hif.Grid.Rows]);
    [~,nm,~]=fileparts(hfl{1});
    nm=regexp(nm,'A(?<year>\d+)(?<day>\d{3})','match');
    nm=cell2mat(nm);
    tfn=sprintf('%s%02i.tif',[wkpth nm '_p'],n);
    matV2tif(tfn,Et,xll,yll,rs_o,ndv,ors_o,wkpth);

% Image processing (read, mosaic, project, crop, resample)
    Tfl=[Tfl;tfn];
  end
  tfl=mat2cell(Tfl,ones(size(Tfl,1),1),size(Tfl,2));
  tfl=regexp(tfl,'(?<year>\d{4})(?<day>\d{3})','match');
  ds=cellfun(@cell2mat,tfl,'UniformOutput',false);
  ds=cell2mat(unique(ds));
  ds=datestr(doy2date(str2double(ds(5:end)),str2double(ds(1:4))),'yyyymmdd');
  a=regexp(vrf,'_\d+m');
  oun=[opth vrf(1:a-1) ds '.tif'];
  imo=MODISimg(Tfl,[],[],wkpth,oun,xl,xr,rx,yb,yt,ry,ors,'bilinear');
  imo(imo==ndv)=NaN;
end
end

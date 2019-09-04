function imout=SRTMimg(fnl,wkpth,nv_o,ndv,thr,ors,xl,xr,rx,yb,yt,ry,oun,itm)
%% Check the input
switch nargin
  case {1:12}; error('Not enough number of arguments');
  case 13; itm='bilinear';
  case 14
  otherwise; error('Too many number of arguments');
end

idn=fullfile(wkpth,['idZ_r' num2str(rx,'%i') '.mat']);
xo=-180;
yo=55;
Trg=5;

for n=1:size(fnl,1)
%% Unzip the file
  fn=unzip(fnl(n,:),wkpth);
  for i=1:length(fn)
    k=strfind(fn{i},'.tif');
    if ~isempty(k)
      fn1=fn{i};
    end
  end

%% Upscale the image
  p=double(imread(fn1));
  p(p==nv_o)=NaN;
  if n==1
    rso=deg2km(Trg/size(p,1),'earth')*1000;
    kx=max(1,fix(rx/rso));
    ky=max(1,fix(ry/rso));
  end
  p=resizeimg(p,ndv,kx,ky,idn,thr);

%% Output the upscaled image
  [~,nm,~]=fileparts(fn1);
  ci=str2double(nm(6:7));
  ri=str2double(nm(9:10));
  xll=xo+Trg*(ci-1);
  yll=yo-Trg*(ri-1);
  tfn=[wkpth 'DEM_' num2str(n,'%02i') '.tif'];
  matV2tif(tfn,p,xll,yll,Trg/size(p,1),ndv,'wgs84',wkpth)
end

%% Build virtual mosaiced image
fun='gdalbuildvrt -overwrite ';
inv=[wkpth 'DEM_*.tif'];
im1=fullfile(wkpth,'DEM.vrt');

system([fun '"' im1 '" ' inv]); % On linux

%% Project image
fun='gdalwarp -of GTiff -overwrite ';
pr1=['-t_srs ' ors ' '];
pr2=sprintf('-te %i %i %i %i ',xl,yb,xr,yt);
pr3=['-r ' itm ' '];
if ~isempty(rx) && ~isempty(ry)
  pr4=sprintf('-tr %i %i ',rx,ry);
else
  pr4=[];
end
prm=[pr1 pr2 pr3 pr4];
ouv=fullfile(wkpth,'DEM_p.tif');

system([fun prm '"' im1 '" "' ouv '"']);

%% Move the processed image to ouv
imout=double(imread(ouv));
if isempty(oun)
  delete(ouv);
else
  movefile(ouv,oun);
end

delete(fullfile(wkpth,'srtm*'));
delete(fullfile(wkpth,'DEM*'));
end

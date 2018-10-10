% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 11/10/2017

%% Functionality
% This function calculates the percentage of a land cover type for a block.

%% Input
%  zo : the input image in high resolution;
%  Xo : the left, right boundary and resolution of X coordinate of input image;
%  Yo : the top, bottum boundary and resolution of Y coordinate of input image;
%  nx : this can be a scale ratio of X coordinate or the left, right boundary
%       and resolution of X coordinate of querry image (for the scale ratio case,
%       if x resolution of the querry image is 10 times coarser, nx is 10);
%  ny : this can be a scale ratio of Y coordinate or the top, bottum boundary
%       and resolution of Y coordinate of querry image;
% infn: full name of the id.mat file (use this to avoid repeatly running knnsearch);
% lcvR: value of a land cover type;
% thr : a threshold percentage of nodata value for a block where block with percent
%       of nodata value greater than this will be assigned the nodata value;
% ndv : a scalar for nodata value.

%% Output
% Zq: the percentage of the land cover type.

function Zq=LC_Frac(zo,Xo,nx,Yo,ny,infn,lcvR,thr,ndv)
if isempty(thr)
  thr=1.0001;
end
zo(isnan(zo))=ndv;

% Search points
if isempty(Xo) && isempty(Yo)
  x=.5:size(zo,2)-.5;
  y=size(zo,1)-.5:-1:.5;
else
  x=Xo(1)+Xo(3)/2:Xo(3):Xo(2)-Xo(3)/2;
  y=Yo(1)-Yo(3)/2:-Yo(3):Yo(2)+Yo(3)/2;
end
[X,Y]=meshgrid(x,y);
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);

% Query points
if isempty(Xo) && isempty(Yo)
  x=nx/2:nx:size(zo,2)-nx/2;
  y=size(zo,1)-ny/2:-ny:ny/2;
else
  if length(nx)==3 && length(ny)==3
    x=nx(1)+nx(3)/2:nx(3):nx(2)-nx(3)/2;
    y=ny(1)-ny(3)/2:-ny(3):ny(2)+ny(3)/2;

    nx=nx(3)/Xo(3);
    ny=ny(3)/Yo(3);

  else
    x=Xo(1)+Xo(3)*nx/2:Xo(3)*nx:Xo(2)-Xo(3)*nx/2;
    y=Yo(1)-Yo(3)*ny/2:-Yo(3)*ny:Yo(2)+Yo(3)*ny/2;
  end
end
[Xq,Yq]=meshgrid(x,y);
Xq=reshape(Xq,numel(Xq),1);
Yq=reshape(Yq,numel(Yq),1);

% Find the N nearest neighbors
N=max([1 round(pi()*nx*ny/4)]); % Search area
if exist(infn,'file')~=0
  load(infn,'id');
else
  [id,~]=knnsearch([X,Y],[Xq,Yq],'k',N);
  save(infn,'id');
end
clear X Y

% Fraction of classes
Z=zo(id);
Z=Z==lcvR;
Zq=sum(Z,2)/N;

Z=zo(id);
Z=Z==ndv;
mk=sum(Z,2)/N; % Mask out location with too many NaN
mk=mk<thr;
Zq(~mk)=ndv;

Zq=reshape(Zq,length(y),length(x));
end

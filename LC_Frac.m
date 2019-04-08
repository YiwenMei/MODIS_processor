% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 08/22/2017

%% Functionality
% This function calculates the percentage of a land cover type for a block.

%% Input
%  Zo : the input image in high resolution;
% ndv : no-data value assigned to the output images;
%  nx : this can be a scale ratio of X coordinate or the left, right boundary
%       and resolution of X coordinate of querry image (for the scale ratio case,
%       if x resolution of the querry image is 10 times coarser, nx is 10);
%  ny : this can be a scale ratio of Y coordinate or the top, bottum boundary
%       and resolution of Y coordinate of querry image;
% idfn: full name of the id.mat file used to avoid repeatly running knnsearch
%       (in line 70) when this function is called wihtin a loop (the code create
%       and save the id.mat in the first iteration and then directly load the
%       id.mat instead of calculating it again to speed up the processing);
% lcvR: value of a land cover type;
% thr : a percentage of nodata value for a block where block with this percent
%       of nodata value is assigned the nodata value;
%  Xo : the left, right boundary and resolution of X coordinate of input image;
%  Yo : the top, bottum boundary and resolution of Y coordinate of input image.

%% Output
% Zq: the percentage of the land cover type.

function Zq=LC_Frac(Zo,ndv,nx,ny,idfn,lcvR,thr,Xo,Yo)
%% Check the input
switch nargin
  case {1:6}; error('Not enough number of arguments');

  case 7
    x=.5:size(Zo,2)-.5;
    y=size(Zo,1)-.5:-1:.5;

    xq=nx/2:nx:size(Zo,2)-nx/2;
    yq=size(Zo,1)-ny/2:-ny:ny/2;

  case 8; error('Y coordinate missing');

  case 9
    x=Xo(1)+Xo(3)/2:Xo(3):Xo(2)-Xo(3)/2;
    y=Yo(1)-Yo(3)/2:-Yo(3):Yo(2)+Yo(3)/2;

    if length(nx)==3 && length(ny)==3
      xq=nx(1)+nx(3)/2:nx(3):nx(2)-nx(3)/2;
      yq=ny(1)-ny(3)/2:-ny(3):ny(2)+ny(3)/2;
      nx=nx(3)/Xo(3);
      ny=ny(3)/Yo(3);

    else
      xq=Xo(1)+Xo(3)*nx/2:Xo(3)*nx:Xo(2)-Xo(3)*nx/2;
      yq=Yo(1)-Yo(3)*ny/2:-Yo(3)*ny:Yo(2)+Yo(3)*ny/2;
    end
    
  otherwise; error('Too many number of arguments');
end

%% Search and query points
[X,Y]=meshgrid(x,y);
X=reshape(X,numel(X),1);
Y=reshape(Y,numel(Y),1);

[Xq,Yq]=meshgrid(xq,yq);
Xq=reshape(Xq,numel(Xq),1);
Yq=reshape(Yq,numel(Yq),1);

%% Find the N nearest neighbors
N=max([1 round(pi()*nx*ny/4)]); % Search area
if exist(idfn,'file')~=0
  load(idfn,'id');
else
  [id,~]=knnsearch([X,Y],[Xq,Yq],'k',N);
  save(idfn,'id');
end
clear X Y

%% Fraction of classes
Z=Zo(id);
Z=Z==lcvR;
Zq=sum(Z,2)/N;

Z=Zo(id);
Z=Z==ndv;
mk=sum(Z,2)/N; % Mask out location with too many NaN
mk=mk<thr;
Zq(~mk)=ndv;

Zq=reshape(Zq,length(yq),length(xq));
end

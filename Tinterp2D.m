function vi=Tinterp2D(vib,via,dn)
[~,dnb,~]=fileparts(vib{1});
dnb=datenum(dnb(3:end),'yyyymmdd'); % datenum of VI record before the inqury time step
vib=read2Dvar(vib);
wb=1-abs((dn-dnb))/8;

[~,dna,~]=fileparts(via{1});
dna=datenum(dna(3:end),'yyyymmdd'); % datenum of VI record after the inqury time step
via=read2Dvar(via);
wa=1-abs((dn-dna))/8;

vi=nansum(cat(3,vib*wb,via*wa),3);

k1=isnan(vib) | isnan(via);
vi(k1)=nansum([vib(k1) via(k1)],2);
vi(isnan(vib) & isnan(via))=NaN; % Final VI with 1-month lag
end

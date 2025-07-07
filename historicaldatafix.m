load('data/updated_historical_isotopes.mat')
%%%d18O - see Hansen et al. 2013
oxyrecord=d18o_y_highfid;
deepoceantemp(oxyrecord<1.75)=(-4*oxyrecord(oxyrecord<1.75))+12;
deepoceantemp(oxyrecord>=1.75)=-2*(oxyrecord(oxyrecord>=1.75)-4.25);
%deepoceantemp=5-(8/3)*(oxyrecord-1.75);
gast_pleistocenemetric=(2*deepoceantemp)+12.25;
gast_pliocenemetric=(2.5*deepoceantemp)+12.15;
anchortime=find(d18o_x_highfid<-5.3,1);
anchorvaldeepocean=deepoceantemp(anchortime);
anchorvalpliocene=gast_pliocenemetric(anchortime);
avgsurftemps_highfid=deepoceantemp-anchorvaldeepocean+anchorvalpliocene;

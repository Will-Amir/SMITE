%%% live plotting script for biomass density

land=double(land_past);
land(land ==0) = NaN ; 

a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(currentcull, [0 20])) ;
a.EdgeColor = "none" ;
clim([0,1]);
title('Biomass growback threshhold');
set(a, 'CData', circshift(currentcull, [0 20]));
colorbar
drawnow
%{
cmap=parula(max(max(biohistories(:,:,end)))+1);
cmap(1,:)=[1,0,0];
%}
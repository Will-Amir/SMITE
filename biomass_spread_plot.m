%%% live plotting script for the spread of plant life

land=double(land_past);
land(land ==0) = NaN ; 


b = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(culledmaterial, [0 20])) ;
b.EdgeColor = "none" ;
%clim([0,0.01]);
title('Biomass spread map');
set(b, 'CData', circshift(culledmaterial, [0 20]));
colorbar
drawnow
%exportgraphics(gca,"grwothspread.gif","Append",true)
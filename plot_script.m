%%plotting script

%download m_proj

land = double(land_past) ; 
land(land ==0) = NaN ; 

figure
subplot(2,3,1)
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift( land, [0 20])) ;
a.EdgeColor = "none" ; 
hold on
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(final_biomass_past, [0 20])) ;
a.EdgeColor = "none" ; 
colorbar
title('Biomass (gC/m^{2})')
hold off

subplot(2,3,2)
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift( land, [0 20])) ;
a.EdgeColor = "none" ; 
hold on
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(NPP_past .* land, [0 20])) ;
a.EdgeColor = "none" ; 
colorbar
title('NPP (gC/m^{2}/year)')
hold off

subplot(2,3,3)
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(double(land_past),[0 20])) ; 
a.EdgeColor = "none" ; 
colorbar
title('Land mask')

subplot(2,3,6)
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift( land, [0 20])) ;
a.EdgeColor = "none" ; 
hold on
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(firemap, [0 20])) ;
a.EdgeColor = "none" ; 
colorbar
title('Wildfire regions')
hold off
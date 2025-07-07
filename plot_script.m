%%plotting script

%download m_proj
global biohistories
land = double(land_past) ; 
land(land ==0) = NaN ; 

figure
subplot(2,3,1)
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift( land, [0 20])) ;
a.EdgeColor = "none" ; 
hold on
cmap=parula(max(max(biohistories(:,:,end)))+1);
cmap(1,:)=[1,0,0];
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(biohistories(:,:,end), [0 20])) ;
colormap(cmap)
a.EdgeColor = "none" ; 
colorbar
%txt = ['Year: ' num2str(t_geol*1e6)];
%text(0,0,txt,'FontSize',16)
title('Biomass (gC/m^{2})')
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

savefig("plantgrowth")
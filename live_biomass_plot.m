%%% live plotting script for biomass density

land=double(land_past);
land(land ==0) = NaN ; 

checkflag=0;
a = pcolor(INTERPSTACK.lon,INTERPSTACK.lat,circshift(currentcull, [0 20])) ;
a.EdgeColor = "none" ;
clim([0,1]);
title('Map of current cullfractions');
set(a, 'CData', circshift(currentcull, [0 20]));
txt = num2str(t_geol*1e6)+"yr";
text(325,0,txt,'FontSize',14,'Color',"white")
colorbar
drawnow
if mod(t_geol*1e6,10)==0
    exportgraphics(gca,"growthspread.gif","Append",true,"Resolution",125)
    checkflag=1;
elseif rem(t_geol*1e6,10)==0
    if checkflag==0
        exportgraphics(gca,"growthspread.gif","Append",true,"Resolution",125)
        checkflag=1;
    end
end
%{
cmap=parula(max(max(biohistories(:,:,end)))+1);
cmap(1,:)=[1,0,0];
%}
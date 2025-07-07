t_geol = 18 ; 
CO2_lvl = 8; 
global CO2ppm
% 1=10 ;     9=360 ;        17=2240 ;       25=56000 ;
% 2=50 ;     10=420 ;       18=2800 ;       26=112000 ;
% 3=100 ;    11=560 ;       19=3360 ; 
% 4=160 ;    12=700 ;       20=4200 ;
% 5 =180 ;   13=840 ;       21=5600 ; 
% 6=200 ;    14=1120 ;      22=7000 ;
% 7=240 ;    15=1400 ;      23=14000 ;
% 8=280 ;    16=11680 ;     24=28000 ;
 
config = 18 ; 
% 1=-540 ;    6=-400 ;       11=-470 ;      16=-145 ;       21 =-15 ;  
% 2=-500 ;    7=-370 ;       12=-245 ;      17=-90 ;        22=0 ; 
% 3=-470 ;    8=-340 ;       13=-220 ;      18=-70 ; 
% 4=-450 ;    9=-300 ;       14=-200 ;      19=-52 ; 
% 5=-430 ;    10=-280 ;      15=-180 ;      20=-30 ; 
 
RUNOFF_past = INTERPSTACK.runoff(:,:,CO2_lvl, config) ;
RUNOFF_future = RUNOFF_past ; 
Tair_past = INTERPSTACK.Tair(:,:,CO2_lvl, config) ; 
Tair_future = Tair_past ; 
land_past = INTERPSTACK.land(:,:,config) ; 
land_future = land_past ; 
TOPO_past = INTERPSTACK.topo(:,:,config) ; 
TOPO_future = TOPO_past ; 
GRID_AREA_km2 = INTERPSTACK.gridarea ; 
pO2 = 20.95 ; 
mrO2 = pO2 / 1000 ; 
contribution_past = 1 ; 
contribution_future = 0 ; 
EVO = 1 ; 
RCO2 = 1 ; 
    
    Flora_roottypes %%%SCRIPT

if runcontrol == -3
    n = 1 ; 
    for i = 1:40
        for j = 1:48
            temp = singlerun.final_biomass(i,j,:) ; 
            grid(:,n+1) = temp(:) ; 
            temp = singlerun.root(i,j,:) ; 
            Croot_grid(:,n+1) = temp(:) ;
            temp = singlerun.stem(i,j,:) ; 
            Cstem_grid(:,n+1) = temp(:) ; 
            temp = singlerun.leaf(i,j,:) ; 
            Cleaf_grid(:,n+1) = temp(:) ; 
            n = n + 1 ; 
        end
    end
    subplot(4,5,7)
    plot(grid)
    title('Final biomass')
    subplot(4,5,8)
    plot(Croot_grid)
    title('Croot')
    subplot(4,5,9)
    plot(Cstem_grid)
    title('C stem')
    subplot(4,5,10)
    plot(Cleaf_grid)
    title('C leaf')

   
    subplot(4,5,11)
    imagesc(singlerun.biome(:,:,end).* land_past)
    colorbar
    title('Biome')
    subplot(4,5,12)
    imagesc(singlerun.ag_biomass(:,:,end) .* land_past)
    colorbar
    title('AG Biomass')
    subplot(4,5,13)
    imagesc(singlerun.bg_biomass(:,:,end) .* land_past )
    colorbar
    title('BG Biomass')

    
    load forcings\INTERPSTACK_2021_improved.mat
    load data/ORNLbiomass.mat
    load data/FLORA_present_day_data.mat
    %all in gC/m^2 
    subplot(4,5,14)
    plot(INTERPSTACK.lat, sum(singlerun.ag_biomass(:,:,end),2,'omitnan'))
    hold on
    plot(lat_cru, sum(flipud(ag_biomass),2,'omitnan') )
    xlabel('Lat')
    title('AG biomass')
    subplot(4,5,15)
    plot(INTERPSTACK.lon, sum(singlerun.ag_biomass(:,:,end),'omitnan'))
    hold on
    plot(lon_cru+ 180, sum(flipud(ag_biomass),'omitnan') )
    xlabel('Lon')
    subplot(4,5,16)
    plot(INTERPSTACK.lat, sum(singlerun.bg_biomass(:,:,end),2,'omitnan'))
    hold on
    plot(lat_cru, sum(flipud(bg_biomass),2,'omitnan') )
    xlabel('Lat')
    title('BG biomass')
    subplot(4,5,17)
    plot(INTERPSTACK.lon, sum(singlerun.bg_biomass(:,:,end),'omitnan'))
    hold on
    plot(lon_cru+180, sum(flipud(bg_biomass),'omitnan') )
    xlabel('Lon')
    legend('Modelled','Actual')

 
    

end


% figure
% subplot(3,1,1)
% imagesc(singlerun.Cleaf)
% title('C_leaf')
% colorbar
% subplot(3,1,2)
% imagesc(singlerun.Croot)
% title('C_root')
% colorbar
% subplot(3,1,3)
% imagesc(singlerun.Cstem)
% title('C_stem')
% colorbar
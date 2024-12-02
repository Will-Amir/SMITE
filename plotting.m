%%%% Proxy color chart
pc1 = [65 195 199]./255 ;
pc2 = [73 167 187]./255 ;
pc3 = [82 144 170]./255 ;
pc4 = [88 119 149]./255 ;
pc5 = [89 96 125]./255 ;
pc6 = [82 56 100]./255 ;
pc7 = [0.9290 0.6940 0.1250] ; 
pc8 = [199 141 247]./ 255 ; 
pc9 = [0.4660 0.6740 0.1880] ;
pc10 = [245 152 231]./255 ; 

load('data\geochem_data_2020.mat')
load('data\Scotese_GAT_2021.mat')
%%% Relative vegetation
figure
hold on
box on
xlabel('Time (Ma)')
ylabel('Relative Vegetation')
load sens_runs.mat
%%%% plot this model
plot(sens.time_myr(351:end,:), nanmean(sens.VEG(351:end,:),2), 'color', c_mean)
plot(sens.time_myr(351:end,:), nanmax(sens.VEG(351:end,:),[],2), 'color',c_range)
plot(sens.time_myr(351:end,:), nanmin(sens.VEG(351:end,:),[],2),'color',c_range)
 
%%% Temp
figure
hold on
plot((sens.time_myr(351:end,:)),nanmean(sens.T_gast(351:end,:),2),'linewidth',1,'color',c_mean)
plot((sens.time_myr(351:end,:)),max(sens.T_gast(351:end,:),[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr(351:end,:)),min(sens.T_gast(351:end,:),[],2),'linewidth',0.5,'color',c_range)
 
%%% CO2 
figure
set(gca, 'YScale', 'log')
hold on
box on
ylim([100 10000])
xlabel('Time (Ma)')
ylabel('Atmospheric CO_{2} (ppm)')
%%%% plot data comparison
%%%% paleosol
% errorbar(paleosol_age,paleosol_co2,paleosol_low,paleosol_high,'color',[0.4 0.7 0.7],'linestyle','none')
plot(paleosol_age(1:534), paleosol_co2(1:534),'s','markeredgecolor',pc1)
%%%% alkenone
% errorbar(alkenone_age,alkenone_co2,alkenone_low,alkenone_high,'color',[0.4 0.7 0.4],'linestyle','none')
plot(alkenone_age, alkenone_co2,'.','markerfacecolor',pc7,'markeredgecolor',pc7)
%%%% boron
% errorbar(boron_age,boron_co2,boron_low,boron_high,'color',[0.7 0.4 0.4],'linestyle','none')
plot(boron_age, boron_co2,'o','markeredgecolor',pc8)
%%%% stomata
% errorbar(stomata_age,stomata_co2,stomata_low,stomata_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(stomata_age(1:271),stomata_co2(1:271),'d','markeredgecolor',pc9)
%%%% liverwort
% errorbar(liverwort_age,liverwort_co2,liverwort_low,liverwort_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(liverwort_age, liverwort_co2,'v','markeredgecolor',pc10, 'markerfacecolor', pc10)
%%%% phytane
% errorbar(phytane_age,phytane_co2,phytane_low,phytane_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(phytane_age(1:255), phytane_co2(1:255),'.','markerfacecolor',pc6,'markeredgecolor',pc6)
%%%% plot this model
plot((sens.time_myr(351:end,:)),nanmean(sens.CO2ppm(351:end,:),2),'linewidth',1,'color',c_mean)
plot((sens.time_myr(351:end,:)),max(sens.CO2ppm(351:end,:),[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr(351:end,:)),min(sens.CO2ppm(351:end,:),[],2),'linewidth',0.5,'color',c_range)
 
%%% Biomass maps + runoff, silicate weathering 
biomass_colour = [ %140 81 10
                   191 127 45
                    223 194 125
                   246 232 195
                   199 234 229
                   128 205 193
                   53 151 143 
                   1 102 94 ] ./ 255 ; 
                   %0 60 48 ] ./ 255 ; 
load('C:\Users\fbskgu\Desktop\Matlab\Weathering\forcings\INTERPSTACK_2021_improved.mat')
lon = INTERPSTACK.lon ; 
lat = INTERPSTACK.lat; 

% timestamps used for Figure 2: 22(0Ma), 17(90Ma), 16(145Ma), 14(200Ma),
% 13(220Ma) 
m_proj( 'robinson', 'longitutde' , [ 0 352.5 ] , 'latitude' , [ -87 87 ] )
land = double( run.gridstate.land(:,:,13) ) ; 
land( land == 0 ) = NaN ;
land( land == 1 ) = -1e5 ; 
m_pcolor( lon, lat, land) 
hold on
m_pcolor( lon, lat, log10( run.gridstate.BIOMASS_tot_grid(:,:,22)))
m_grid 
colormap( biomass_colour )
caxis([2.5 4.25])
colorbar
hold off
title( 'Ordovician: 470 Ma' )
 
runoff_colour = [209 229 240
                146 197 222
                107 175 214
                49 130 189
                8 81 156] ./ 240 ; 
land = double( run.gridstate.land(:,:,13) ) ; 
land( land == 0 ) = NaN ;
m_pcolor( lon, lat, run.gridstate.Q(:,:,13) .* land)
m_grid
colormap(runoff_colour)
colorbar
caxis([0 1000])
 
sil_colour = [%222 233 250
            189 210 242
            153 156 207
            145 136 189
            119 93 153] ./ 255 ;
 
m_pcolor(lon,lat,run.gridstate.CW(:,:,13) .* run.gridstate.f_biota(:,:,13))
colormap(sil_colour)
colorbar
caxis([0 12])
m_grid('fontsize', 0.1)
 

for i = 1:22
    habitable_area(i) = sum(sum(~isnan(run.gridstate.BIOMASS_tot_grid(:,:,i)) .* area)) ;
    total_area(i) = sum(nansum(run.gridstate.land(:,:,i) .* area)) ;
    rel_area(i) = habitable_area(i)/total_area(i) *100 ;
    sil_weathering(i) = sum(sum(~isnan(run.gridstate.CW(:,:,i)) .* area)) ;
    weathering(i) = sum(sum(~isnan(run.gridstate.f_biota(:,:,i)) .* area)) ;
end
time = run.gridstate.time_myr(11:end) ; 
figure
plot(time, weathering(11:end)/weathering(end))
hold on
plot(time, sil_weathering(11:end)/sil_weathering(end))
ylabel('Relative weathering')
xlabel('Time (Ma)')
legend('fbiota', 'silicate', 'Location', 'southeast')

figure
plot(rel_area(11:end))
ylabel('Habitable area (%)')
xlabel('Time (Ma)')
%%%minimum threshold - 632 gC/m2


%%%% d13C record
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot data comparison
plot(d13c_x,d13c_y,'.','color',pc2)
%%%% plot this model
plot(run.state.time_myr,run.state.delta_mccb,'k')
 
%%%% d34S record
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('\delta^{34}S_{sw}')
%%%% plot data comparison
plot(d34s_x,d34s_y,'.','color',pc2)
%%%% plot this model
plot(run.state.time_myr,run.state.d34s_S,'k')

%%%% Ocean 87Sr/86Sr 
hold on
box on
xlim([-250 run.pars.whenend/1e6])
ylim([0.706 0.71])
xlabel('Time (Ma)')
ylabel('^{87}Sr/^{86}Sr seawater')
%%%% plot data comparison
plot(sr_x,sr_y,'color',pc2)
%%%% plot this model
plot(run.state.time_myr,run.state.delta_OSr,'k')

%%%% SO4
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Marine SO_{4} (mM)')
%%%% plot algeo data window comparison
plot(sconc_max_x,sconc_max_y,'color',pc1)
plot(sconc_min_x,sconc_min_y,'color',pc1)
plot(sconc_mid_x,sconc_mid_y,'color',pc2)
%%%% plot fluid inclusion data comparison
for u = 1:2:length(SO4_x-1)
   plot( [SO4_x(u) SO4_x(u)] , [SO4_y(u) SO4_y(u+1)], 'color' , pc3 ) ;     
end
%%%% plot this model
plot(run.state.time_myr,(run.state.S./run.pars.S0)*28,'k')

%%%% O2 (%) 
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Atmospheric O_{2} (%)')
%%%% plot data comparison
for u = 1:2:length(O2_x) - 1
   plot( [O2_x(u) O2_x(u)] , [O2_y(u) O2_y(u+1)] , 'color' , pc2  ) ;     
end
%%%% plot this model
plot(run.state.time_myr,run.state.mrO2.*100,'k')

%%%% ICE LINE
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Ice line')
%%%% plot iceline proxy
plot(paleolat_x,paleolat_y,'color' ,pc1) ;
%%%% plot this model
plot(run.state.time_myr,run.state.iceline,'k') ;
ylim([40 90])
colormap(gca,'gray')

%%% Corg fluxes
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Flux (mol/yr)')
%%%% plot this model
plot(run.state.time_myr,run.state.mocb,'b')
plot(run.state.time_myr,run.state.locb,'g')
plot(run.state.time_myr,run.state.oxidw,'r')
plot(run.state.time_myr,run.state.ocdeg,'k') 
%%%% Legend
text(-250,5e12,'mocb','color','b')
text(-250,4e12,'locb','color','g')
text(-250,3e12,'oxidw','color','r')
text(-250,2e12,'ocdeg','color','k')
%%%% Title
title('C_{org} fluxes')

%%% Ccarb fluxes
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Flux (mol/yr)')
%%%% plot this model
plot(run.state.time_myr,run.state.silw,'r')
plot(run.state.time_myr,run.state.carbw,'c')
plot(run.state.time_myr,run.state.sfw,'b')
plot(run.state.time_myr,run.state.mccb,'k') 
%%%% Legend
text(-250,28e12,'silw','color','r')
text(-250,24e12,'carbw','color','c')
text(-250,20e12,'sfw','color','b')
text(-250,16e12,'mccb','color','k')
%%%% Title
title('C_{carb} fluxes')

%%% S fluxes
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
% ylim([0 5e12])
ylabel('Fluxes (mol/yr)')
%%%% plot this model
plot(run.state.time_myr,run.state.mpsb,'k')
plot(run.state.time_myr,run.state.mgsb,'c')
plot(run.state.time_myr,run.state.pyrw,'r')
plot(run.state.time_myr,run.state.pyrdeg,'m') 
plot(run.state.time_myr,run.state.gypw,'b')
plot(run.state.time_myr,run.state.gypdeg,'g') 
%%%% Legend
text(-250,1.9e11,'mpsb','color','k')
text(-250,1.7e11,'mgsb','color','c')
text(-250,1.5e12,'pyrw','color','r')
text(-250,1.2e12,'pyrdeg','color','m')
text(-250,1e12,'gypw','color','b')
text(-250,0.8e12,'gypdeg','color','g')
%%%% Title
title('S fluxes')

%%% NUTRIENTS P N
hold on
box on
xlim([-250 run.pars.whenend/1e6])
ylim([0 3])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(run.state.time_myr,run.state.P/run.pars.P0,'b')
plot(run.state.time_myr,run.state.N/run.pars.N0,'g')
%%%% Legend
text(-590,1.5,'P','color','b')
text(-590,1,'N','color','g')
%%%% Title
title('Nutrient reservoirs')

%%%% C SPECIES
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(run.state.time_myr,run.state.G/run.pars.G0,'k')
plot(run.state.time_myr,run.state.C/run.pars.C0,'c')
plot(run.state.time_myr,run.state.VEG,'g--')
%%%% Legend
text(-590,1.5,'VEG','color','g')
text(-590,1.25,'G','color','k')
text(-590,1,'C','color','c')
%%%% Title
title('C reservoirs')

%%%% S SPECIES
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(run.state.time_myr,run.state.PYR/run.pars.PYR0,'k')
plot(run.state.time_myr,run.state.GYP/run.pars.GYP0,'c')
%%%% Legend
text(-590,1,'PYR','color','k')
text(-590,0.9,'GYP','color','c')
%%%% Title
title('S reservoirs')

%%%% Forg and Fpy ratos
hold on
box on
xlim([-250 run.pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('f_{org}, f_{py}')
%%%% plot this model
plot(run.state.time_myr,run.state.mocb ./ (run.state.mocb + run.state.mccb),'k')
%%%% plot fpy
plot(run.state.time_myr, run.state.mpsb ./ (run.state.mpsb + run.state.mgsb),'m')

%%%% GLOBAL FORCINGS
hold on
box on
xlim([-250 run.pars.whenend/1e6])
ylim([0 2.5])
xlabel('Time (Ma)')
ylabel('Relative forcing')
%%%% plot this model
plot(run.state.time_myr,run.state.DEGASS,'r')
plot(run.state.time_myr,run.state.BAS_AREA,'k')
plot(run.state.time_myr,run.state.EVO,'g')
plot(run.state.time_myr,run.state.W,'b')
plot(run.state.time_myr,run.state.Bforcing,'m')
plot(run.state.time_myr,run.state.GRAN_AREA,'color',[0.8 0.8 0.8])
%%%% Legend
text(-590,2.4,'D','color','r')
text(-590,2.2,'E','color','g')
text(-590,2,'W','color','b')
text(-590,1.8,'B','color','m')
text(-590,1.6,'BA','color','k')
text(-590,1.4,'GA','color',[0.8 0.8 0.8])
%%%% Title
title('Forcings')

%%% Colour scheme for Figure S2
biomass_colour = [ 255 255 204
    220 232 149
    161 218 180
    132 219 161
    65 182 196
    33 149 163
    44 127 184
    22 104 161
    37 52 148
    11 20 79 ] ./ 255 ; 

grayscale = flipud( [ 20 20 20
    37 38 38
    59 57 57
    92 88 88
    107 103 103
    125 120 120
    143 137 137
    171 164 164
    194 188 188
    212 206 205 
    230 223 223
    230 230 230 ]./ 255) ;

%Calculating difference in data curves: 250Ma onwards; binned every Ma
%modelled data
temp = scifi_run.state.tempC(889:end) ; 
time = scifi_run.state.time_myr(889:end) ;
%SCION data
scion_time = scion_run.state.time_myr(600:end) ; 
scion_temp = scion_run.state.tempC(600:end) ; 
%actual data interpolated to fit model time
time_bins = -250:0 ; 
interp_temp(:,1) = interp1(flipud(Scotese_2021_age(1:251)), flipud(Scotese_2021_GAT(1:251)), time_bins) ;
%SCION data interpolated to fit model time
interp_temp(:,2) = interp1(scion_time, scion_temp, time_bins) ; 
interp_temp(:,3) = interp1(time, temp, time_bins) ; 
data_diff = interp_temp(:,1) - interp_temp(:,2) ;
data_diff(:,2) = interp_temp(:,1) - interp_temp(:,3) ; 
dist = sqrt(data_diff(:,1).^2);
dist(:,2) = sqrt(data_diff(:,2).^2) ; 


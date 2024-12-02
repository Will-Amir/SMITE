%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    SCI-FI: SCION-FLORA Biomass comparison to CEDA       %%%
%%%                 Present day validation                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load data/FLORA_present_day_data.mat land_cru lat_cru lon_cru runoff tmp_avg
load data/area_present.mat
 %%%%%% loading parameters for biomass calculations
 
% Longitude and latitude 
x_lon = 360 ; 
y_lat = 720 ; 

% Photosynthesis parameters (taken from LPJ model)
t25 = 2600 ; 
q_10t = 0.57 ; 
v = 0.7 ; 
kc25 = 30 ; 
ko25 = 3e4 ; 
q_10c = 2.1 ; 
q_10o = 1.2 ; 
s = (24/12) * 0.015 ; %0.08 in Github for C3 plants
alpha = 0.08 ; 
theta = 0.7 ; 
lr_max = 0.75 ; 
CN_leaf = 29 ;
% Tissue respiration rate at 10 degree C
r_tem = 0.055 * 365 ; % gC/gN/d -> gC/gN/year
r_bor = 0.066 * 365 ;
r_tro = 0.011 * 365 ; 
% Growth respiraition
R_growth = 0.25 ;     
% Leaf longevity; ranges from 0.5 - 1 depending on type of plant
life_leaf_tem = 0.75 ; 
life_leaf_bor = 0.75 ;
life_leaf_tro = 1 ; 
% Minimum weathering
minbiota = 0.32 ; 

%%%%%%% Insolation
ins = 150 + 250 .* normpdf(lat_cru, 0 , 40 ) ./ normpdf( 0, 0, 40 ) .* ones( x_lon, y_lat ) ;

%%%%%%% Ice (< -10 degC) or no runoff = no biomass 
tmp_avg( tmp_avg < -20 ) = NaN ; 

water_stress = runoff ; 
water_stress(water_stress == 0 | isnan(water_stress)) = NaN ; 
water_stress = 1 - (1 ./ ( 1 + exp(0.005 .* (water_stress - 450)))) ;
pO2 = 20.9 * 1000 ; 
%%%%%%% Photosynthesis calculation for the past keyframe
tf = t25 * ( q_10t .^ ( ( tmp_avg - 25 ) * 0.1 ) ) ; 
tstar = pO2 ./ ( 2 * tf ) ; 

intra_pp = v * 280 ; 

kc = kc25 * ( q_10c .^ ( ( tmp_avg - 25 ) * 0.1 ) ) ; 
ko = ko25 * ( q_10o .^ ( ( tmp_avg - 25 ) * 0.1 ) ) ; 


c2 = ( intra_pp - tstar ) ./ ( intra_pp + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 

ftemp_tem = normpdf( tmp_avg, 15,7 ) ; % 15, 7 ) ; % Temperate 
ftemp_bor = normpdf( tmp_avg, 0, 20 ) ; %5, 10 ) ; % Boreal 
ftemp_tro = normpdf( tmp_avg, 27, 7 ) ; % 27, 3 ) ;  % Tropical 

c1_tem = alpha * ftemp_tem .* ( ( intra_pp - tstar ) ./ ( intra_pp + 2 * tstar ) ) ;
c1_bor = alpha * ftemp_bor .* ( ( intra_pp - tstar ) ./ ( intra_pp + 2 * tstar ) ) ; 
c1_tro = alpha * ftemp_tro .* ( ( intra_pp - tstar ) ./ ( intra_pp + 2 * tstar ) ) ;

sigma = ( 1 - ( ( c2 - s ) ./ ( c2 - theta * s ) ) ) .^ 0.5 ; 

g_T = exp( 308.56 .* ( ( 1 / 56.02 ) - ( 1 ./ ( tmp_avg + 46.02 ) ) ) ) ;

photosynth_tem = 10 .* 365 * ins .* ( c1_tem ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* water_stress ; 
photosynth_bor = 10 .* 365 * ins .* ( c1_bor ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* water_stress ; 
photosynth_tro = 10 .* 365 * ins .* ( c1_tro ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* water_stress ; 

%%%%%%% Carbon in leaf allocation - first timestep
C_leaf_tem = lr_max .* photosynth_tem ; 
C_leaf_bor = lr_max .* photosynth_bor ;
C_leaf_tro = lr_max .* photosynth_tro ; 


%%%%%%% Biomass starting at homogenous values
biomass_tem( :, :, 1 ) = 2.5e4 * land_cru ; 
biomass_tem( biomass_tem == 0 ) = NaN ; 
biomass_bor( :, :, 1 ) = 2.5e4 * land_cru ; 
biomass_bor( biomass_bor == 0 ) = NaN ; 
biomass_tro( :, :, 1 ) = 2.5e4 * land_cru ; 
biomass_tro( biomass_tro == 0 ) = NaN ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Calculating Biomass   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Stopping at <1% change in biomass
biomass_change_final_step(1) = 999 ; 
n = 1 ; 
biomass_tot(1) = 0 ; 
%%%%%%% Turnover changes between 8% and 20% depending on O2
turnover = min( max( 0.0092 *( pO2 - 10 ), 0.08 ), 0.2 ) ; 

while abs( biomass_change_final_step ) > 1 

    %%% Leaf respiration per plant type
    R_leaf_tem = r_tem * ( C_leaf_tem / CN_leaf ) .* g_T ;
    R_leaf_bor = r_bor * ( C_leaf_bor / CN_leaf ) .* g_T ;
    R_leaf_tro = r_tro * ( C_leaf_tro / CN_leaf ) .* g_T ;
    R_leaf_tem( R_leaf_tem <= 0 ) = 0 ; 
    R_leaf_bor( R_leaf_bor <= 0 ) = 0 ; 
    R_leaf_tro( R_leaf_tro <= 0 ) = 0 ; 

    %%% NPP
    NPP_tem = ( 1 - R_growth ) .* ( photosynth_tem - R_leaf_tem ) ; 
    NPP_bor = ( 1 - R_growth ) .* ( photosynth_bor - R_leaf_bor ) ; 
    NPP_tro = ( 1 - R_growth ) .* ( photosynth_tro - R_leaf_tro ) ; 

    %%% Biomass
    biomass_tem = biomass_tem + ( C_leaf_tem - turnover * biomass_tem )  ; 
    biomass_bor = biomass_bor + ( C_leaf_bor - turnover * biomass_bor )  ; 
    biomass_tro = biomass_tro + ( C_leaf_tro - turnover * biomass_tro )  ; 
    final_biomass = max( biomass_tem, max( biomass_bor, biomass_tro ) )  ;

    %%% Carbon in leaf allocation (n+1)
    C_leaf_tem = ( C_leaf_tem .* ( 1 - life_leaf_tem ) ) + ( lr_max .* NPP_tem ) ;
    C_leaf_bor = ( C_leaf_bor .* ( 1 - life_leaf_bor ) ) + ( lr_max .* NPP_bor ) ;
    C_leaf_tro = ( C_leaf_tro .* ( 1 - life_leaf_tro ) ) + ( lr_max .* NPP_tro ) ;

    %%% Biomass total (gC/m2)
    biomass_tot(n+1) = sum( nansum( final_biomass .* area_present ) ) ; 
    biomass_change_final_step = ( ( biomass_tot( n + 1 ) - biomass_tot( n ) ) / biomass_tot( n + 1 ) ) * 100 ; 

    n = n + 1 ; 
end

%%%%%%%%%%%%%%%%%%%% END OF RUN %%%%%%%%%%%%%%%%%%%%

%%%% CEDA Biomass data (2020) %%%%
%Original tiff files taken from: https://catalogue.ceda.ac.uk/uuid/bedc59f37c9545c981a839eb552e4084
%29 tiff files; original resolution = 45000x45000
%Scaled original files to 45x45 

load data\CEDA_biomass.mat
lat_ceda(:,1) = [-89:90] ;
lon_ceda(:,1) = [-180:0.89:180] ;
%Combining 29 subsets into one map (9 columns; 4 rows)
%Blank = missing gridcells for the oceans

% CEDA_biomass = [n80w180 n80w140 n80w100 n80w060 n80w020 n80e020 n80e060 n80e100 n80e140 ; 
%                 n40w180 n40w140 n40w100 n40w060 n40w020 n40e020 n40e060 n40e100 n40e140 ; 
%                 n00w180 blank   n00w100 n00w060 n00w020 n00e020 n00e060 n00e100 n00e140 ; 
%                 blank   blank   s40w100 s40w060 blank   blank   blank   blank   s40e140 ] ; 

%Dataset resolution: 180 x 405
%Model resolution: 40 x 48
%Comparing model to actual biomass
%CEDA biomass = Mg/ha
%Model biomass = g/m^2
%Mg/ha * 100 -> g/m^2 **Both maps are now in g/m^2**
% CEDA_low_res = imresize(CEDA_biomass * 100 , [40 48]) ; 
CEDA_low_res = double(CEDA_low_res) ; 
load data\model_present.mat
load forcings\INTERPSTACK_2021_improved.mat 
lon = INTERPSTACK.lon ; 
lat = INTERPSTACK.lat ; 
land(land == 0) = NaN ; 
CEDA_low_res = CEDA_low_res .* land ; 
compare_colour = [ 84 48 5 
                   140 81 10 
                   191 129 45
                   223 194 125
                   246 232 195
                   245 245 245
                   199 234 229 
                   128 205 193 
                   53 151 143 
                   1 102 94 
                   0 60 48 ] ./ 246 ;

biomass_colour = [  223 194 125
                   140 81 10
                   191 127 45
                   246 232 195
                   199 234 229
                   128 205 193
                   53 151 143 
                   1 102 94
                   0 60 48 ] ./ 246 ; 

%Requires m_proj mapping package
m_proj( 'robinson', 'longitutde' , [ 0 352.5 ] , 'latitude' , [ -87 87 ] )
figure
subplot(2,1,1)
hold on
m_pcolor(lon,lat,land)
CEDA_low_res(CEDA_low_res== 0 ) = NaN ; 
m_pcolor(lon,lat,CEDA_low_res)
m_grid
colormap(biomass_colour)
a = colorbar ; 
ylabel(a,'gC/m^{2}','FontSize',10,'Rotation',360);
caxis([0 1e4])
title('CEDA biomass')
hold off

subplot(2,1,2)
hold on
m_pcolor(lon,lat,land)
m_pcolor(lon,lat,model_present)
m_grid
a = colorbar ; 
ylabel(a,'gC/m^{2}','FontSize',10,'Rotation',360);
caxis([0 1e4])
title('Model biomass')
hold off

figure
hold on
m_pcolor(lon,lat,land)
m_pcolor(lon,lat,CEDA_low_res - model_present)
m_grid
colormap(compare_colour)
a = colorbar ; 
ylabel(a,'gC/m^{2}','FontSize',10,'Rotation',360);
caxis([-1e4 1e4])
title('Difference in biomass')

figure
subplot(2,1,1)
plot( lon, nansum( CEDA_low_res ) )
hold on
plot( lon, nansum( model_present ) )
ylabel( 'g C/ m^{2}' )
title ( 'Longitude' )

subplot(2,1,2)
plot( lat, nansum( CEDA_low_res , 2 ) )
hold on
plot( lat, nansum( model_present , 2 ) ) 
ylabel( 'g C/ m^{2}' ) 
title( 'Latitude' ) 
legend( 'Actual data', 'Model data' )
hold off

figure
m_pcolor(lon,lat,land)
hold on
model_ice = double(model_temp <-10) ; 
model_ice(model_ice ==1 ) = -10 ; 
model_ice(model_ice == 0 ) = NaN ; 
m_pcolor(lon,lat,model_ice)
m_grid
title('Ice covered land')
function dy = SCIFI_equations2(t,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FLORA - Fast Land Occupancy Reaction Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BJW Mills 2021 (SCION)
%%%% b.mills@leeds.ac.uk
%%%% KG 2023 (FLORA)
%%%% k.gurung@leeds.ac.uk

%%%% model equations

%%%%%%% setup dy array
dy = zeros(21,1);  
%%%%%%% set up global parameters
global stepnumber
global pars
global forcings
global workingstate
global gridstate
global INTERPSTACK
global sensanal
global sensparams
global biopars
global random_impactor_flag
global asteroidtimes
global lats
global longs
global powers
global numasteroids
global timetoinject
global culledmaterial
global SMITEflag
global seeding_primer

%%%%%%% get variables from Y to make working easier
P = y(1) ;
O = y(2) ;
A = y(3) ;
S = y(4) ;
G = y(5) ;
C = y(6) ;
PYR = y(7) ;
GYP = y(8) ;
N = y(11) ;
OSr = y(18) ;
SSr = y(20) ;
dSSr = y(21)/y(20) ;

%%%%%%% geological time in Ma
t_geol = t*(1e-6) ;

%%%%%%% calculate isotopic fractionation of reservoirs
delta_G = y(12)/y(5);
delta_C = y(13)/y(6);
delta_GYP  = y(15)/y(8);
delta_PYR  = y(14)/y(7);

%%%%%%% atmospheric fraction of total CO2, atfrac(A)
atfrac0 = 0.01614 ; % constant
atfrac = atfrac0 * (A/pars.A0) ; % variable

%%%%%%% calculations for pCO2, pO2
RCO2 = (A/pars.A0)*(atfrac/atfrac0) ;
CO2atm = RCO2*(280e-6) ;
CO2ppm = RCO2*280 ;

%%%%%%% mixing ratio of oxygen (not proportional to O reservoir)
mrO2 = ( O/pars.O0 )  /   ( (O/pars.O0)  + pars.copsek16 ) ;
%%%%%%% relative moles of oxygen 
RO2 =  O/pars.O0 ;
%%%%%%% pO2 = mixing ratio * atmospheric pressure
pO2 = mrO2 * 1000 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   Interpolate forcings for this timestep   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% COPSE Reloaded forcing set
E_reloaded = interp1( 1e6 * forcings.t, forcings.E , t ) ;
W_reloaded = interp1( 1e6 * forcings.t, forcings.W , t ) ;
%%%%%%% Additional forcings
GR_BA = interp1( forcings.GR_BA(:,1) , forcings.GR_BA(:,2) , t ) ;
newGA = interp1( forcings.newGA(:,1) , forcings.newGA(:,2) , t ) ;
D_combined_mid = interp1( forcings.D_force_x , forcings.D_force_mid, t_geol) ;
D_combined_min = interp1( forcings.D_force_x , forcings.D_force_min, t_geol) ;
D_combined_max = interp1( forcings.D_force_x , forcings.D_force_max, t_geol) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Choose forcing functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEGASS = D_combined_mid ;
W = W_reloaded ;
EVO = E_reloaded ;
CPLAND = 1 ;
Bforcing = interp1qr([-1000 -150 -100 0]',[0.75 0.75 1 1]',t_geol) ;
BAS_AREA = GR_BA ;
GRAN_AREA = newGA ;
capdelS = 27   ;
capdelC_land = 27   ;
capdelC_marine = 35  ;

%%%%%%% SHORELINE
SHORELINE = interp1qr(forcings.shoreline_time',forcings.shoreline_relative',t_geol) ;

%%%%%%% bioturbation forcing
f_biot = interp1qr([-1000e6 -525e6 -520e6 0]',[0 0 1 1]',t);
CB = interp1qr([0 1]',[1.2 1]',f_biot) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Sensitivity analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% all sensparams vary between [-1 +1]
if sensanal == 1
    
    %%% Very degassing between upper and lower bounds
    if sensparams.randminusplus1 > 0
        DEGASS = (1 - sensparams.randminusplus1)*DEGASS + sensparams.randminusplus1*D_combined_max ;
    else
        DEGASS = (1 + sensparams.randminusplus1)*DEGASS - sensparams.randminusplus1*D_combined_min ;
    end
        
    %%% simple +/- 20% variation
    BAS_AREA = BAS_AREA * (1 + 0.2*sensparams.randminusplus2) ;
    GRAN_AREA = GRAN_AREA * (1 + 0.2*sensparams.randminusplus3) ;
   
    %%%
    capdelS = 30 + 10*sensparams.randminusplus5 ;
    capdelC_land = 25 + 5*sensparams.randminusplus6 ;
    capdelC_marine = 30 + 5*sensparams.randminusplus7 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial fields from stack   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Fetch keyframe grids   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% find past and future keyframes
key_past_time = max ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) <= 0 ) ) ;
key_future_time = min ( INTERPSTACK.time( (INTERPSTACK.time - t_geol) >= 0 ) ) ;
if isempty(key_past_time) == 1
    key_past_time = key_future_time ;
end

%%%%%%% find keyframe indexes and fractional contribution
key_past_index = find( INTERPSTACK.time == key_past_time ) ;
key_future_index = find( INTERPSTACK.time == key_future_time ) ;
dist_to_past = abs( key_past_time - t_geol ) ;
dist_to_future = abs( key_future_time - t_geol ) ;
%%%%%%% fractional contribution of each keyframe
if dist_to_past + dist_to_future == 0
    contribution_past = 1 ;
    contribution_future = 0 ;
else
    contribution_past = dist_to_future / ( dist_to_past + dist_to_future ) ;
    contribution_future = dist_to_past / ( dist_to_past + dist_to_future ) ;
end

%%%%%%% intrepolate keyframe CO2 concentrations to generate keyframe fields
%%%%%%% find INTERPSTACK keyframes using model CO2
key_upper_CO2 = min ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) >= 0 ) ) ;
key_lower_CO2 = max ( INTERPSTACK.CO2( (INTERPSTACK.CO2 - CO2ppm) <= 0 ) ) ;
%%%%%%% find keyframe indexes and fractional contribution
key_upper_CO2_index = find( INTERPSTACK.CO2 == key_upper_CO2 )  ;
key_lower_CO2_index = find( INTERPSTACK.CO2 == key_lower_CO2 ) ;
dist_to_upper = abs( key_upper_CO2 - CO2ppm ) ;
dist_to_lower = abs( key_lower_CO2 - CO2ppm ) ;

%%%%%%% fractional contribution of each keyframe
if dist_to_upper + dist_to_lower == 0
    contribution_lower = 1 ;
    contribution_upper = 0 ;
else
    contribution_upper = dist_to_lower / ( dist_to_upper + dist_to_lower ) ;
    contribution_lower = dist_to_upper / ( dist_to_upper + dist_to_lower ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create time keyframes using CO2 keyfield contributions   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Runoff
RUNOFF_past = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_past_index); 
RUNOFF_future = contribution_upper.*INTERPSTACK.runoff(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.runoff(:,:,key_lower_CO2_index,key_future_index); 
%%%%%%% Tair
Tair_past = contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_past_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_past_index); 
Tair_future = contribution_upper.*INTERPSTACK.Tair(:,:,key_upper_CO2_index,key_future_index) + contribution_lower.*INTERPSTACK.Tair(:,:,key_lower_CO2_index,key_future_index); 

%%%%%%% time kayframes that don't depend on CO2
%%%%%%% Topography
TOPO_past = INTERPSTACK.topo(:,:,key_past_index) ; 
TOPO_future = INTERPSTACK.topo(:,:,key_future_index) ; 

%%%%%%% last keyframe land recorded for plot
land_past = INTERPSTACK.land(:,:,key_past_index) ; 
land_future = INTERPSTACK.land(:,:,key_future_index) ; 

%%%%%%% gridbox area
GRID_AREA_km2 = INTERPSTACK.gridarea ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Global variables   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% global average surface temperature
GAST = mean(mean( Tair_past .* pars.rel_contrib ))*contribution_past  +  mean(mean( Tair_future .* pars.rel_contrib ))*contribution_future  ;

%%%%%%% basalt and granite temp dependency - direct and runoff
Tsurf = GAST + 273 ;
TEMP_gast = Tsurf ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  SMITE - Sudden Magnitudinous Impactor-Triggered Evolution  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
struck = INTERPSTACK.struck(:,:) ;
potentialfires = INTERPSTACK.wildfires(:,:) ;
burntmaterial = zeros(40,48) ;
if random_impactor_flag==0
    if any(asteroidtimes(:)<=t_geol)
        for time = 1:length(asteroidtimes)
            if asteroidtimes(time)<=t_geol
                valstorun=[asteroidtimes(time),lats(time),longs(time),powers(time)];
                SMITE
                SMITEflag=1;
                asteroidtimes(time)=1;
            end
        end
   end
else
    if t_geol>=timetoinject
        SMITE
        SMITEflag=1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   FLORA - Fast Land Occupant Reaction Alogrithm   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Insolation
ins = biopars.ins_present - ( biopars.ins_present * 4.6/100 * ( abs( t_geol )/570 ) ) ; 

%%%%%%% Ice (< -10 degC) or no runoff = no biomass 
Tair_ice_past = Tair_past ; 
Tair_ice_future = Tair_future ; 
Tair_ice_past( Tair_ice_past < -30 ) = NaN ; 
Tair_ice_future( Tair_ice_future < -30 ) = NaN ; 

% water_stress_past = zeros( biopars.x_lon, biopars.y_lat ) ; 
% water_stress_future = zeros( biopars.x_lon, biopars.y_lat ) ; 
water_stress_past = RUNOFF_past ; 
water_stress_past(water_stress_past == 0 | isnan(water_stress_past)) = NaN ; 
water_stress_future = RUNOFF_future ; 
water_stress_future(water_stress_future == 0 | isnan(water_stress_future)) = NaN ; 
water_stress_past = 1 - (1 ./ ( 1 + exp(0.005 .* (water_stress_past - 450)))) ;
water_stress_future = 1 - (1 ./ ( 1 + exp(0.005 .* (water_stress_future - 450)))) ;
% for i = 1 : biopars.x_lon
%     for j = 1 : biopars.y_lat
%         if RUNOFF_past( i, j ) == 0 
%             water_stress_past( i, j ) = NaN ; 
%         else
%             water_stress_past( i, j ) = sigmf( RUNOFF_past( i, j ) , [0.005 450] ) ; 
%         end
% 
%         if RUNOFF_future( i, j ) == 0 
%             water_stress_future( i, j ) = NaN ; 
%         else
%             water_stress_future( i, j ) = sigmf( RUNOFF_future( i, j ) , [0.005 450] ) ; 
%         end
%     end
% end

%%%%%%% Photosynthesis calculation for the past keyframe
tf = biopars.t25 * ( biopars.q_10t .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 
tstar = pO2 ./ ( 2 * tf ) ; 

pi = biopars.v * CO2ppm ; 

kc = biopars.kc25 * ( biopars.q_10c .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 
ko = biopars.ko25 * ( biopars.q_10o .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 

c2 = ( pi - tstar ) ./ ( pi + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 

ftemp_tem = normpdf( Tair_ice_past, 15,7 ) ; % 15, 7 ) ; % Temperate 
ftemp_bor = normpdf( Tair_ice_past, 0, 20 ) ; %5, 10 ) ; % Boreal 
ftemp_tro = normpdf( Tair_ice_past, 27, 7 ) ; % 27, 3 ) ;  % Tropical 

c1_tem = biopars.alpha * ftemp_tem .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;
c1_bor = biopars.alpha * ftemp_bor .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ; 
c1_tro = biopars.alpha * ftemp_tro .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;

sigma = ( 1 - ( ( c2 - biopars.s ) ./ ( c2 - biopars.theta * biopars.s ) ) ) .^ 0.5 ; 

g_T_past = exp( 308.56 .* ( ( 1 / 56.02 ) - ( 1 ./ ( Tair_ice_past + 46.02 ) ) ) ) ;

photosynth_tem_past = 10 .* 365 * ins .* ( c1_tem ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_past ; 
photosynth_bor_past = 10 .* 365 * ins .* ( c1_bor ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_past ; 
photosynth_tro_past = 10 .* 365 * ins .* ( c1_tro ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_past ; 

%%%%%%% Carbon in leaf allocation - first timestep
C_leaf_tem_past = biopars.lr_max .* photosynth_tem_past ; 
C_leaf_bor_past = biopars.lr_max .* photosynth_bor_past ;
C_leaf_tro_past = biopars.lr_max .* photosynth_tro_past ; 

%%%%%%% Photosynthesis calculation for the future keyframe
tf = biopars.t25 * ( biopars.q_10t .^ ( ( Tair_ice_future - 25 ) * 0.1 ) ) ; 
tstar = pO2 ./ ( 2 * tf ) ; 

pi = biopars.v * key_lower_CO2 ; 

kc = biopars.kc25 * ( biopars.q_10c .^ ( ( Tair_ice_future - 25 ) * 0.1 ) ) ; 
ko = biopars.ko25 * ( biopars.q_10o .^ ( ( Tair_ice_future - 25 ) * 0.1 ) ) ; 

c2 = ( pi - tstar ) ./ ( pi + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 

ftemp_tem = normpdf( Tair_ice_future, 15,7 ) ; %15, 7 ) ; % Temperate 
ftemp_bor = normpdf( Tair_ice_future, 0, 20 ) ; %5, 10 ) ; % Boreal 
ftemp_tro = normpdf( Tair_ice_future, 27, 7 ) ; %27, 3 ) ;  % Tropical 

c1_tem = biopars.alpha * ftemp_tem .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;
c1_bor = biopars.alpha * ftemp_bor .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ; 
c1_tro = biopars.alpha * ftemp_tro .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;

sigma = ( 1 - ( ( c2 - biopars.s ) ./ ( c2 - biopars.theta * biopars.s ) ) ) .^ 0.5 ; 

g_T_future = exp( 308.56 .* ( ( 1 / 56.02 ) - ( 1 ./ ( Tair_ice_future + 46.02 ) ) ) ) ;

photosynth_tem_future = 10 .* 365 * ins .* ( c1_tem ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_future ; 
photosynth_bor_future = 10 .* 365 * ins .* ( c1_bor ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_future ; 
photosynth_tro_future = 10 .* 365 * ins .* ( c1_tro ./ c2 ) .* ( c2 - ( ( ( 2 * biopars.theta ) -1 ) * biopars.s ) - ( 2 .* ( c2 - ( biopars.theta .* biopars.s ) ) .* sigma ) ) .* water_stress_future ; 

%%%%%%% Carbon in leaf allocation - first timestep
C_leaf_tem_future( :, :, 1 ) = biopars.lr_max .* photosynth_tem_future ; 
C_leaf_bor_future( :, :, 1 ) = biopars.lr_max .* photosynth_bor_future ;
C_leaf_tro_future( :, :, 1 ) = biopars.lr_max .* photosynth_tro_future ; 

%%%%%%% Seeding initial biomass starting at homogenous values
if seeding_primer == 0
    biomass_tem_past=workingstate.biomass;
    biomass_bor_past=workingstate.biomass;
    biomass_tro_past=workingstate.biomass;
    biomass_tem_future=workingstate.biomass;
    biomass_bor_future=workingstate.biomass;
    biomass_tro_future=workingstate.biomass;
else
    biomass_tem_past( :, :, 1 ) = 2.5e4*land_past;%.*culledmaterial ; 
    biomass_tem_past( biomass_tem_past == 0 ) = NaN ; 
    biomass_bor_past( :, :, 1 ) = 2.5e4*land_past;%.*culledmaterial ; 
    biomass_bor_past( biomass_bor_past == 0 ) = NaN ; 
    biomass_tro_past( :, :, 1 ) = 2.5e4*land_past;%.*culledmaterial ; 
    biomass_tro_past( biomass_tro_past == 0 ) = NaN ;  
    biomass_tem_future( :, :, 1 ) = 2.5e4*land_future;%.*culledmaterial ; 
    biomass_tem_future( biomass_tem_future == 0 ) = NaN ; 
    biomass_bor_future( :, :, 1 ) = 2.5e4*land_future;%.*culledmaterial ; 
    biomass_bor_future( biomass_bor_future == 0 ) = NaN ; 
    biomass_tro_future( :, :, 1 ) = 2.5e4*land_future;%.*culledmaterial ; 
    biomass_tro_future( biomass_tro_future == 0 ) = NaN ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%   Calculating Biomass   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Stopping at <1% change in biomass
biomass_change_final_step(1) = 999 ; 
n = 1 ; 

%%%%%%% Turnover changes between 8% and 20% depending on O2
turnover = min( max( 0.0092 *( mrO2 * 100 - 10 ), 0.08 ), 0.2 ) ; 

while abs( biomass_change_final_step ) > 1 

    %%% Leaf respiration per plant type
    R_leaf_tem_past = biopars.r_tem * ( C_leaf_tem_past / biopars.CN_leaf ) .* g_T_past ;
    R_leaf_bor_past = biopars.r_bor * ( C_leaf_bor_past / biopars.CN_leaf ) .* g_T_past ;
    R_leaf_tro_past = biopars.r_tro * ( C_leaf_tro_past / biopars.CN_leaf ) .* g_T_past ;
    R_leaf_tem_future = biopars.r_tem * ( C_leaf_tem_future / biopars.CN_leaf ) .* g_T_future ;
    R_leaf_bor_future = biopars.r_bor * ( C_leaf_bor_future / biopars.CN_leaf ) .* g_T_future ;
    R_leaf_tro_future = biopars.r_tro * ( C_leaf_tro_future / biopars.CN_leaf ) .* g_T_future ;
    R_leaf_tem_past( R_leaf_tem_past <= 0 ) = 0 ; 
    R_leaf_bor_past( R_leaf_bor_past <= 0 ) = 0 ; 
    R_leaf_tro_past( R_leaf_tro_past <= 0 ) = 0 ; 
    R_leaf_tem_future( R_leaf_tem_future <= 0 ) = 0 ; 
    R_leaf_bor_future( R_leaf_bor_future <= 0 ) = 0 ; 
    R_leaf_tro_future( R_leaf_tro_future <= 0 ) = 0 ; 

    %%% NPP
    NPP_tem_past = ( 1 - biopars.R_growth ) .* ( photosynth_tem_past - R_leaf_tem_past ) ; 
    NPP_bor_past = ( 1 - biopars.R_growth ) .* ( photosynth_bor_past - R_leaf_bor_past ) ; 
    NPP_tro_past = ( 1 - biopars.R_growth ) .* ( photosynth_tro_past - R_leaf_tro_past ) ; 
    NPP_tem_future = ( 1 - biopars.R_growth ) .* ( photosynth_tem_future - R_leaf_tem_future ) ;
    NPP_bor_future = ( 1 - biopars.R_growth ) .* ( photosynth_bor_future - R_leaf_bor_future ) ;
    NPP_tro_future = ( 1 - biopars.R_growth ) .* ( photosynth_tro_future - R_leaf_tro_future ) ;

    %%% Biomass
    biomass_tem_past = biomass_tem_past + ( C_leaf_tem_past - turnover * biomass_tem_past) * biopars.dt ; 
    biomass_bor_past = biomass_bor_past + ( C_leaf_bor_past - turnover * biomass_bor_past ) * biopars.dt ; 
    biomass_tro_past = biomass_tro_past+ ( C_leaf_tro_past - turnover * biomass_tro_past ) * biopars.dt ; 
    biomass_tem_future = biomass_tem_future + ( C_leaf_tem_future - turnover * biomass_tem_future ) * biopars.dt ;
    biomass_bor_future = biomass_bor_future + ( C_leaf_bor_future - turnover * biomass_bor_future ) * biopars.dt ;
    biomass_tro_future = biomass_tro_future + ( C_leaf_tro_future - turnover * biomass_tro_future ) * biopars.dt ;

    final_biomass_past = max( biomass_tem_past, max( biomass_bor_past, biomass_tro_past ) ) ;
    final_biomass_future = max( biomass_tem_future, max( biomass_bor_future, biomass_tro_future ) ) ; 
%     NPP_past = max( NPP_tem_past, max( NPP_bor_past, NPP_tro_past ) ) * EVO ; 
%     NPP_future = max( NPP_tem_future, max( NPP_bor_future, NPP_tro_future ) ) * EVO ; 

    %%% Carbon in leaf allocation (n+1)
    C_leaf_tem_past = ( C_leaf_tem_past .* ( 1 - biopars.life_leaf_tem ) ) + ( biopars.lr_max .* NPP_tem_past ) ;
    C_leaf_bor_past = ( C_leaf_bor_past .* ( 1 - biopars.life_leaf_bor ) ) + ( biopars.lr_max .* NPP_bor_past ) ;
    C_leaf_tro_past = ( C_leaf_tro_past .* ( 1 - biopars.life_leaf_tro ) ) + ( biopars.lr_max .* NPP_tro_past ) ;
    C_leaf_tem_future = ( C_leaf_tem_future .* ( 1 - biopars.life_leaf_tem ) ) + ( biopars.lr_max .* NPP_tem_future ) ;
    C_leaf_bor_future = ( C_leaf_bor_future .* ( 1 - biopars.life_leaf_bor ) ) + ( biopars.lr_max .* NPP_bor_future ) ;
    C_leaf_tro_future = ( C_leaf_tro_future .* ( 1 - biopars.life_leaf_tro ) ) + ( biopars.lr_max .* NPP_tro_future ) ;


    %%% Biomass total (gC/m2)
    biomass_past_tot = sum( sum( final_biomass_past .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' ) ) ; 
    biomass_future_tot = sum( sum( final_biomass_future .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' ) ) ; 
    biomass_tot( n + 1 ) = biomass_past_tot * contribution_past + biomass_future_tot * contribution_future ; 
    biomass_change_final_step = ( ( biomass_tot( n + 1 ) - biomass_tot( n ) ) / biomass_tot( n + 1 ) ) * 100 ; 

    n = n + 1 ; 
end

%%%%%%% Calculating effect of biomass on weathering
for i = 1 : 40
    for j = 1 : 48
        %{
       %%%%%%% Mild random biomass reduction effects
       if SMITEflag==1
           toadd=abs((t_geol-timetoinject))/20;
           culledmaterial(i,j)=min(1,(culledmaterial(i,j)+toadd));
       end
        %}
       %1 = temperate, 2 = boreal, 3 = tropical, 4 = ice/desert
       if final_biomass_past ( i , j ) == biomass_tem_past( i , j ) * EVO
            biome( i , j ) = 1 ; 
           NPP_past(i,j) = NPP_tem_past(i,j);%*culledmaterial(i,j) ; 
       elseif final_biomass_past(i,j) == biomass_bor_past( i , j ) * EVO
            biome( i , j ) = 2 ; 
           NPP_past(i,j) = NPP_bor_past(i,j);%*culledmaterial(i,j) ; 
       elseif final_biomass_past(i,j) == biomass_tro_past( i , j ) * EVO
            biome( i , j ) = 3 ;
           NPP_past(i,j) = NPP_tro_past(i,j);%*culledmaterial(i,j) ; 
       else
            biome( i , j ) = 4 ; 
           NPP_past(i,j) = NaN ; 
       end
       if final_biomass_future ( i, j ) == biomass_tem_future( i, j ) * EVO
           NPP_future( i, j ) = NPP_tem_future( i, j );
       elseif final_biomass_future ( i, j ) == biomass_bor_future( i, j ) * EVO
           NPP_future( i, j ) = NPP_bor_future( i, j );
       elseif final_biomass_future ( i, j ) == biomass_tro_future( i, j ) * EVO
           NPP_future( i, j ) = NPP_tro_future( i, j );
       else
           NPP_future( i, j ) = NaN ; 
       end
       %%%%%%% Calculate the culling effects of an asteroid
       final_biomass_past(i,j)=(final_biomass_past(i,j));%*(culledmaterial(i,j)));
       final_biomass_future(i,j)=(final_biomass_future(i,j));
       %burntmaterial(i,j)=(final_biomass_past(i,j)*(1-culledmaterial(i,j)))/100;
    end
end

f_biota_past = 0.0005 * NPP_past * EVO + ( biopars.minbiota * RCO2 ^ 0.25 ) ; 
f_biota_future = 0.0005 * NPP_future * EVO + ( biopars.minbiota * RCO2 ^ 0.25 ) ; 

f_biota_past( land_past == 1 & isnan( f_biota_past ) ) = biopars.minbiota .* RCO2 ^ 0.25  ;
f_biota_future( land_future == 1 & isnan( f_biota_future ) ) = biopars.minbiota .* RCO2 ^ 0.25 ; 
workingstate.biomass=final_biomass_past;
%disp(workingstate.biomass)
%%%%%%% Mass of biosphere
VEG = biomass_tot( end ) / 4.53e17 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   SMITE CONTINUED   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firecarbon=(sum(sum( sum( burntmaterial .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' ) ))/12);
%*interp1([-70 -66.2 -66.1 -66 -65.9 -65.8 -60],[0 0 1 1 1 0 0],t_geol);
dl13c_firecarbon=-5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Spatial silicate weathering   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Equation from Joshi et al., 2019 GRL

%%%%%%% runoff in mm/yr
Q_past = RUNOFF_past ;
Q_past(Q_past<0) = 0 ;
Q_future = RUNOFF_future ;
Q_future(Q_future<0) = 0 ;

%%%%%%% temp in kelvin
T_past = Tair_past + 273 ;
T_future = Tair_future + 273 ;

%%%%%%% Erosion rate from Maffre et al. 2018
model_elevation_past = TOPO_past ;
model_elevation_future = TOPO_future ;
model_elevation_past(isnan(model_elevation_past)) = 0 ;
model_elevation_future(isnan(model_elevation_future)) = 0 ;

%%%%%%% make lat and lon grids for erosion calculation
lat_grid = INTERPSTACK.lat'.*ones(40,48) ;
lon_grid = INTERPSTACK.lon.*ones(40,48) ;
%%%%%%% gradient calculation
[~,~,dFdyNorth_past,dFdxEast_past] = gradientm(lat_grid,lon_grid,model_elevation_past) ;
[~,~,dFdyNorth_future,dFdxEast_future] = gradientm(lat_grid,lon_grid,model_elevation_future) ;
%%%%%%% topographic slope
tslope_past = ( dFdyNorth_past.^2 + dFdxEast_past.^2 ).^0.5 ;
tslope_future = ( dFdyNorth_future.^2 + dFdxEast_future.^2 ).^0.5 ;
%%%%%%% pierre erosion calculation, t/m2/yr
k_erosion = 3.3e-3 ; %%%% for 16Gt present day erosion in FOAM
EPSILON_past = k_erosion .* (Q_past.^0.31) .* tslope_past .* max(Tair_past,2) ;
EPSILON_future = k_erosion .* (Q_future.^0.31) .* tslope_future .* max(Tair_future,2) ;
%%%%%%% check total tonnes of erosion - should be ~16Gt
EPSILON_per_gridbox_past = EPSILON_past .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
EPSILON_per_gridbox_future = EPSILON_future .* GRID_AREA_km2 .* 1e6 ; %%% t/m2/yr * m2
erosion_tot_past = sum(sum(EPSILON_per_gridbox_past)) ;
erosion_tot_future = sum(sum(EPSILON_per_gridbox_future)) ;
erosion_tot = erosion_tot_past*contribution_past + erosion_tot_future*contribution_future ;

%%%%%%% Pierre weathering equation params
Xm = 0.1 ;
K = 6e-5 ; 
kw = 1e-3 ;
Ea = 20 ; 
z = 10 ; 
sigplus1 = 0.9 ; 
T0 = 286 ;
R = 8.31e-3 ;

%%%%%%% equations
R_T_past = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_past) ) ) ;
R_T_future = exp( ( Ea ./ (R.*T0) ) - ( Ea ./ (R.*T_future) ) ) ;
R_Q_past = 1 - exp( -1.*kw .* Q_past ) ;
R_Q_future = 1 - exp( -1.*kw .* Q_future ) ;
R_reg_past = ( (z./EPSILON_past).^sigplus1 ) ./ sigplus1 ;
R_reg_future = ( (z./EPSILON_future).^sigplus1 ) ./ sigplus1 ;

%%%%%%% equation for CW per km2 in each box; grid
CW_per_km2_past = 1e6 .* EPSILON_past .* Xm .* ( 1 - exp( -1.* K .* R_Q_past .* R_T_past .* R_reg_past ) ) ; 
CW_per_km2_future = 1e6 .* EPSILON_future .* Xm .* ( 1 - exp( -1.* K .* R_Q_future .* R_T_future .* R_reg_future ) ) ; 
%%%%%%% CW total
CW_past = CW_per_km2_past .* GRID_AREA_km2 .* f_biota_past ; % f_biota grid 
CW_future = CW_per_km2_future .* GRID_AREA_km2 .* f_biota_future  ;
%%%%%%% world CW
CW_sum_past = sum( CW_past( ~isnan( CW_past ) ) ) ;
CW_sum_future = sum( CW_future( ~isnan( CW_future ) ) ) ;
CW_tot = CW_sum_past*contribution_past + CW_sum_future*contribution_future ;

%%%%%%% carbonate weathering spatial approximation, linear with runoff
k_carb_scale = 200 ; % scaling parameter to recover preent day rate
CWcarb_per_km2_past = k_carb_scale * Q_past ; 
CWcarb_per_km2_future = k_carb_scale * Q_future ; 
%%%%%%% CW total
CWcarb_past = CWcarb_per_km2_past .* GRID_AREA_km2 .* f_biota_past ;
CWcarb_future = CWcarb_per_km2_future .* GRID_AREA_km2 .* f_biota_future ;
%%%%%%% world CWcarb
CWcarb_sum_past = sum( CWcarb_past( ~isnan( CWcarb_past ) ) ) ;
CWcarb_sum_future = sum( CWcarb_future( ~isnan( CWcarb_future ) ) ) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Grid interpolated variables   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% silicate weathering scale factor by present day rate in cation tonnes
silw_scale = 4.2e8 ; %%%% for k erosion 3.3e-3
%%%%%%% overall spatial weathering
silw_spatial = CW_tot * ( (pars.k_basw + pars.k_granw) / silw_scale) ;
carbw_spatial = ( CWcarb_sum_past*contribution_past + CWcarb_sum_future*contribution_future ) ;

%%%%%%% set assumed ice temperature
Tcrit = -10 ;
%%%%%%% ice line calculations
Tair_past_ice = Tair_past ;
Tair_past_ice(Tair_past_ice >= Tcrit) = 0 ;
Tair_past_ice(Tair_past_ice < Tcrit) = 1 ;
Tair_future_ice = Tair_future ;
Tair_future_ice(Tair_future_ice >= Tcrit) = 0 ;
Tair_future_ice(Tair_future_ice < Tcrit) = 1 ;
%%%%%%% count only continental ice
Tair_past_ice = Tair_past_ice.* land_past ;
Tair_future_ice = Tair_future_ice.* land_future ; 
%%%%%%% sum into lat bands
latbands_past = sum(Tair_past_ice,2) ;
latbands_future = sum(Tair_future_ice,2) ;
latbands_past(latbands_past>0) = 1 ;
latbands_future(latbands_future>0) = 1 ;
%%%%%%% find appropiate lat
latresults_past =  INTERPSTACK.lat .* latbands_past' ;
latresults_future =  INTERPSTACK.lat .* latbands_future' ;
latresults_past(latresults_past == 0) = 90 ;
latresults_future(latresults_future == 0) = 90 ;
%%%%%%% lowest glacial latitude
iceline_past = min(abs(latresults_past)) ;
iceline_future = min(abs(latresults_future)) ;
iceline = iceline_past * contribution_past +  iceline_future * contribution_future ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% version using gran area and conserving total silw
basw = silw_spatial * ( pars.basfrac * BAS_AREA / ( pars.basfrac * BAS_AREA + ( 1 - pars.basfrac  ) *GRAN_AREA ) ) ;
granw = silw_spatial * ( ( 1 - pars.basfrac  ) *GRAN_AREA / ( pars.basfrac * BAS_AREA + ( 1 - pars.basfrac  ) *GRAN_AREA ) ) ;
carbw = carbw_spatial ;

%%%%%%% overall weathering
silw = basw + granw ;
carbw_relative = (carbw/pars.k_carbw) ;

%%%%%%% oxidative weathering 
oxidw = pars.k_oxidw*carbw_relative*(G/pars.G0)*((O/pars.O0)^pars.a) ;

%%%%%%% pyrite weathering 
pyrw = pars.k_pyrw*carbw_relative*(PYR/pars.PYR0)  ;

%%%%%%% gypsum weathering 
gypw = pars.k_gypw*(GYP/pars.GYP0)*carbw_relative ;

%%%%%%% seafloor weathering, revised following Brady and Gislason but not directly linking to CO2whenstart
f_T_sfw = exp(0.0608*(Tsurf-288)) ; 
sfw = pars.k_sfw * f_T_sfw * DEGASS ; %%% assume spreading rate follows degassing here

%%%%%%% Degassing 
ocdeg = pars.k_ocdeg*DEGASS*(G/pars.G0) ;
ccdeg = pars.k_ccdeg*DEGASS*(C/pars.C0)*Bforcing ;
pyrdeg = pars.k_pyrdeg*(PYR/pars.PYR0)*DEGASS;
gypdeg = pars.k_gypdeg*(GYP/pars.GYP0)*DEGASS;

%%%%%%% gypsum burial
% mgsb = pars.k_mgsb*(S/pars.S0);
mgsb = pars.k_mgsb*(S/pars.S0)*(1/SHORELINE) ;

%%%%%%% carbonate burial
mccb = carbw + silw ;

%%%%%%% COPSE reloaded P weathering
pfrac_silw = 0.8 ;
pfrac_carbw = 0.14 ;
pfrac_oxidw = 0.06 ;
phosw = pars.k_phosw*( (pfrac_silw)*( silw/pars.k_silw )  +   (pfrac_carbw)*( carbw/pars.k_carbw ) +  (pfrac_oxidw)*(  oxidw/ pars.k_oxidw )  )  ;

%%%%%%% COPSE reloaded
pland = pars.k_landfrac * VEG * phosw  ;
pland0 = pars.k_landfrac*pars.k_phosw;
psea = phosw - pland ;

%%%%%%% convert total reservoir moles to micromoles/kg concentration    
Pconc = ( P/pars.P0 ) * 2.2 ;
Nconc = ( N/pars.N0 ) * 30.9 ;
newp = 117 * min(Nconc/16,Pconc) ;    

%%%%%%% carbon burial
mocb = pars.k_mocb*((newp/pars.newp0)^pars.b) * CB ;
locb = pars.k_locb*(pland/pland0)*CPLAND  ;

%%%%%%% PYR burial function (COPSE)
fox= 1/(O/pars.O0) ; 
%%%%%%% mpsb scales with mocb so no extra uplift dependence
mpsb = pars.k_mpsb*(S/pars.S0)*fox*(mocb/pars.k_mocb) ;

%%%%%%% OCEAN ANOXIC FRACTION
k_anox = 12 ; 
k_u = 0.5 ;
ANOX = 1 / ( 1 + exp( -1 * k_anox * ( k_u * (newp/pars.newp0) - (O/pars.O0) ) ) ) ;

%%%%%%% nutrient burial
CNsea = 37.5 ;
monb = mocb/CNsea ;

%%%%%%% P burial with bioturbation on
CPbiot = 250 ;
CPlam = 1000 ;
mopb = mocb*( (f_biot/CPbiot) + ( (1-f_biot)/CPlam ) ) ;
capb = pars.k_capb*( mocb/pars.k_mocb ) ;

%%%%%%% reloaded
fepb = (pars.k_fepb/pars.k_oxfrac)*(1-ANOX)*(P/pars.P0) ;

%%%%%%% nitrogen cycle
%%%%%%% COPSE reloaded
if (N/16) < P
    nfix = pars.k_nfix * ( ( ( P - (N/16)  ) / (  pars.P0 - (pars.N0/16)    ) )^2 ) ;
else
    nfix = 0 ;
end

denit = pars.k_denit * ( 1 + ( ANOX / (1-pars.k_oxfrac) )  ) * (N/pars.N0) ;

%%%%%%% reductant input
reductant_input = pars.k_reductant_input * DEGASS ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Reservoir calculations  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Phosphate
dy(1) = psea - mopb - capb - fepb ;

%%% Oxygen
dy(2) = locb + mocb - oxidw  - ocdeg  + 2*(mpsb - pyrw  - pyrdeg) - reductant_input ;

%%% Carbon dioxide
dy(3) = -locb - mocb + oxidw + ocdeg + ccdeg + carbw - mccb - sfw  + reductant_input + firecarbon ;

%%% Sulphate
dy(4) = gypw + pyrw - mgsb - mpsb + gypdeg + pyrdeg ;

%%%Buried organic C
dy(5) = locb + mocb - oxidw - ocdeg ;

%%% Buried carb C 
dy(6) = mccb + sfw - carbw - ccdeg ;

%%% Buried pyrite S
dy(7) = mpsb - pyrw - pyrdeg ;

%%% Buried gypsum S 
dy(8) = mgsb - gypw - gypdeg ;

%%%% Nitrate
dy(11) = nfix - denit - monb;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Isotope reservoirs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% d13c and d34s for forwards model
d13c_A = y(16) / y(3) ;
d34s_S = y(17) / y(4) ;

%%%% carbonate fractionation
delta_locb = d13c_A - capdelC_land ; 
delta_mocb = d13c_A - capdelC_marine ; 
delta_mccb = d13c_A ;

%%%%% S isotopes (copse)
delta_mpsb = d34s_S - capdelS ;

%%% deltaORG_C*ORG_C 
dy(12) =  locb*(  delta_locb ) + mocb*( delta_mocb )  -   oxidw*delta_G  -   ocdeg*delta_G  ;

%%% deltaCARB_C*CARB_C 
dy(13) =  mccb*delta_mccb + sfw*delta_mccb  -  carbw*delta_C  - ccdeg*delta_C  ;

%%% deltaPYR_S*PYR_S (young)
dy(14) =  mpsb*( delta_mpsb )  - pyrw*delta_PYR  - pyrdeg*delta_PYR ;

%%% deltaGYP_S*GYP_S (young)
dy(15) =  mgsb*d34s_S   - gypw*delta_GYP  - gypdeg*delta_GYP ;

%%% delta_A * A
dy(16) = -locb*(  delta_locb ) -mocb*( delta_mocb ) + oxidw*delta_G + ocdeg*delta_G + ccdeg*delta_C + carbw*delta_C - mccb*delta_mccb - sfw*delta_mccb + reductant_input*-5 + firecarbon*dl13c_firecarbon;% +noise*dl13c_firecarbon;

%%% delta_S * S
dy(17) = gypw*delta_GYP + pyrw*delta_PYR -mgsb*d34s_S - mpsb*( delta_mpsb ) + gypdeg*delta_GYP + pyrdeg*delta_PYR ; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Strontium system   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% fluxes
Sr_granw = pars.k_Sr_granw *( granw / pars.k_granw ) ;
Sr_basw = pars.k_Sr_basw *( basw / pars.k_basw ) ;
Sr_sedw = pars.k_Sr_sedw *( carbw / pars.k_carbw ) * (SSr/pars.SSr0) ;
Sr_mantle = pars.k_Sr_mantle * DEGASS ;
Sr_sfw = pars.k_Sr_sfw * (sfw/pars.k_sfw) * ( OSr/pars.OSr0 ) ;
Sr_metam = pars.k_Sr_metam * DEGASS * (SSr/pars.SSr0) ;
Sr_sedb = pars.k_Sr_sedb * ( mccb/pars.k_mccb ) * ( OSr/pars.OSr0 ) ;

%%%%%%% fractionation calculations
delta_OSr = y(19) / y(18) ;
delta_SSr = y(21) / y(20) ;

%%%%%%% original frac
RbSr_bas = 0.1 ;
RbSr_gran = 0.26 ;
RbSr_mantle = 0.066 ;
RbSr_carbonate = 0.5 ;

%%%%%%% frac calcs
dSr0 = 0.69898 ;
tforwards = 4.5e9 + t ;
lambda = 1.4e-11 ;
dSr_bas = dSr0 + RbSr_bas*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_gran = dSr0 + RbSr_gran*( 1 - exp(-1*lambda*tforwards) ) ;
dSr_mantle = dSr0 + RbSr_mantle*( 1 - exp(-1*lambda*tforwards) ) ;

%%%%%%% Ocean [Sr]
dy(18) = Sr_granw + Sr_basw + Sr_sedw + Sr_mantle - Sr_sedb - Sr_sfw ;

%%%%%%% Ocean [Sr]*87/86Sr
dy(19) = Sr_granw*dSr_gran + Sr_basw*dSr_bas + Sr_sedw*delta_SSr + Sr_mantle*dSr_mantle - Sr_sedb*delta_OSr - Sr_sfw*delta_OSr ;

%%%%%%% Sediment [Sr]
dy(20) = Sr_sedb - Sr_sedw - Sr_metam ;

%%%%%%% Sediment [Sr]*87/86Sr
dy(21) = Sr_sedb*delta_OSr - Sr_sedw*delta_SSr - Sr_metam*delta_SSr + SSr*lambda*RbSr_carbonate*exp(lambda*tforwards)  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Mass conservation check   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_C = A + G + C ;
res_S = S + PYR + GYP ;
iso_res_C = A*d13c_A + G*delta_G + C*delta_C ;
iso_res_S = S*d34s_S + PYR*delta_PYR + GYP*delta_GYP ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Print full states for single run   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensanal == 0
    workingstate.iso_res_C(stepnumber,1) = iso_res_C ;
    workingstate.iso_res_S(stepnumber,1) = iso_res_S ;
    workingstate.res_C(stepnumber,1) = res_C ;
    workingstate.res_S(stepnumber,1) = res_S ;
    workingstate.time(stepnumber,1) = t;
    workingstate.temperature(stepnumber,1) = TEMP_gast ;
    workingstate.tempC(stepnumber,1) = TEMP_gast - 273 ;
    workingstate.P(stepnumber,1) = P ;
    workingstate.O(stepnumber,1) = O ;
    workingstate.A(stepnumber,1) = A ;
    workingstate.S(stepnumber,1) = S ;
    workingstate.G(stepnumber,1) = G ;
    workingstate.C(stepnumber,1) = C ;
    workingstate.PYR(stepnumber,1) = PYR ;
    workingstate.GYP(stepnumber,1) = GYP ;
    workingstate.N(stepnumber,1) = N ;
    workingstate.OSr(stepnumber,1) = OSr ;
    workingstate.SSr(stepnumber,1) = SSr ;
    %%%%%%% print isotope information
    workingstate.d13c_A(stepnumber,1) = d13c_A ;
    workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
    workingstate.d34s_S(stepnumber,1) = d34s_S ;
    workingstate.delta_G(stepnumber,1) = delta_G ;
    workingstate.delta_C(stepnumber,1) = delta_C ;
    workingstate.delta_PYR(stepnumber,1) = delta_PYR ;
    workingstate.delta_GYP(stepnumber,1) = delta_GYP ;
    workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
    %%%%%%% print forcings
    workingstate.DEGASS(stepnumber,1) = DEGASS ;
    workingstate.W(stepnumber,1) = W ;
    workingstate.EVO(stepnumber,1) = EVO ;
    workingstate.CPLAND(stepnumber,1) = CPLAND ;
    workingstate.Bforcing(stepnumber,1) = Bforcing ;
    workingstate.BAS_AREA(stepnumber,1) = BAS_AREA ;
    workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA ;
    %%%%%%%% print variables
    workingstate.RCO2(stepnumber,1) = RCO2 ;
    workingstate.RO2(stepnumber,1) = RO2 ;
    workingstate.mrO2(stepnumber,1) = mrO2 ;
    workingstate.VEG(stepnumber,1) = VEG ;
    workingstate.ANOX(stepnumber,1) = ANOX ;
    workingstate.iceline(stepnumber,1) = iceline ;
    %%%%%%%% print fluxes
    workingstate.mocb(stepnumber,1) = mocb ;
    workingstate.locb(stepnumber,1) = locb ;
    workingstate.mccb(stepnumber,1) = mccb ;
    workingstate.mpsb(stepnumber,1) = mpsb ;
    workingstate.mgsb(stepnumber,1) = mgsb ;
    workingstate.silw(stepnumber,1) = silw ;
    workingstate.carbw(stepnumber,1) = carbw ;
    workingstate.oxidw(stepnumber,1) = oxidw ;
    workingstate.basw(stepnumber,1) = basw ;
    workingstate.granw(stepnumber,1) = granw ;
    workingstate.phosw(stepnumber,1) = phosw ;
    workingstate.psea(stepnumber,1) = psea ;
    workingstate.nfix(stepnumber,1) = nfix ;
    workingstate.denit(stepnumber,1) = denit ;
    workingstate.pyrw(stepnumber,1) = pyrw ;
    workingstate.gypw(stepnumber,1) = gypw ;
    workingstate.ocdeg(stepnumber,1) = ocdeg ;
    workingstate.ccdeg(stepnumber,1) = ccdeg ;
    workingstate.pyrdeg(stepnumber,1) = pyrdeg ;
    workingstate.gypdeg(stepnumber,1) = gypdeg ;
    workingstate.sfw(stepnumber,1) = sfw ;
    workingstate.Sr_granw(stepnumber,1) = Sr_granw ;
    workingstate.Sr_basw(stepnumber,1) = Sr_basw ;
    workingstate.Sr_sedw(stepnumber,1) = Sr_sedw ;
    workingstate.Sr_mantle(stepnumber,1) = Sr_mantle ;
    workingstate.dSSr(stepnumber,1) = dSSr ;
    workingstate.relativenewp(stepnumber,1) = newp/pars.newp0 ;
    workingstate.erosion_tot(stepnumber,1) = erosion_tot ;
    workingstate.biomass_tot(stepnumber,1) = sum( sum( biomass_tot( end ), 'omitnan' ) ) ;
    workingstate.NPP(stepnumber,1) = ( sum( sum( NPP_past .* ( GRID_AREA_km2 .* 1e6 ), 'omitnan' ) ) * contribution_past ) + ( sum( sum( NPP_future .* (GRID_AREA_km2 .* 1e6 ), 'omitnan' ) ) * contribution_future) ; 
    workingstate.hab_area(stepnumber,1) = ((( sum( sum( ~isnan( final_biomass_past ) & final_biomass_past > 0 ))) / ( sum( sum( land_past, 'omitnan' )))) * 100 * contribution_past ) + ((( sum( sum( ~isnan( final_biomass_future ) & final_biomass_future > 0 )))/ ( sum( sum(land_future, 'omitnan')))) * 100 * contribution_future ) ;
    workingstate.hab_area_raw(stepnumber,1) = sum(sum(~isnan(final_biomass_past) .*GRID_AREA_km2 .*1e6 .*contribution_past, 'omitnan')) + sum(sum(~isnan(final_biomass_future) .* GRID_AREA_km2 .* 1e6 .* contribution_future, 'omitnan')) ; 
%     workingstate.turnover(stepnumber,1) = turnover ; 
    %%%%%%%% print time
    workingstate.time_myr(stepnumber,1) = t_geol ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Print gridstates for single run   %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sensanal == 0
    %%%% print a gridstate when each keytime threshold is crossed, or at model end
    next_stamp = pars.next_gridstamp ;
    if pars.finishgrid == 0
        if t_geol > next_stamp || t_geol == 0

            %%%% write gridstates
            gridstate.time_myr(pars.gridstamp_number,1) = next_stamp ;
            gridstate.land(:,:,pars.gridstamp_number) = land_past ;
            gridstate.Q(:,:,pars.gridstamp_number) = Q_past ;
            gridstate.Tair(:,:,pars.gridstamp_number) = Tair_past ;
            gridstate.TOPO(:,:,pars.gridstamp_number) = TOPO_past ;
            gridstate.CW(:,:,pars.gridstamp_number) = CW_per_km2_past ; %%% t/km2/yr
            gridstate.CWcarb(:,:,pars.gridstamp_number) = CWcarb_past ;
            gridstate.EPSILON(:,:,pars.gridstamp_number) = EPSILON_past * 1e6 ; %%% t/km2/yr
            gridstate.BIOMASS_tot_grid(:,:,pars.gridstamp_number) = final_biomass_past; %g C/m2 
            gridstate.f_biota(:,:,pars.gridstamp_number) = f_biota_past ;
            gridstate.rel_fbiota(:,:,pars.gridstamp_number) = f_biota_past ./ (biopars.minbiota * RCO2 ^ 0.25) ; 
%             gridstate.NPP(:,:,pars.gridstamp_number) = NPP_past ; 
          

             gridstate.biome(:,:,pars.gridstamp_number) = biome ; 
%             gridstate.NPP_biome(:,:,pars.gridstamp_number) = final_NPP ; 
            %%%% set next boundary
            if t_geol < 0
                pars.gridstamp_number = pars.gridstamp_number + 1 ;
                pars.next_gridstamp = pars.runstamps(pars.gridstamp_number) ;
            else
                pars.finishgrid = 1 ;
            end

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Print plotting states only in sensanal   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sensanal == 1
        workingstate.BAS_AREA(stepnumber,1) = BAS_AREA;
        workingstate.GRAN_AREA(stepnumber,1) = GRAN_AREA;
        workingstate.DEGASS(stepnumber,1) = DEGASS;
        workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
        workingstate.d34s_S(stepnumber,1) = d34s_S ;
        workingstate.delta_OSr(stepnumber,1) = delta_OSr ;
        workingstate.SmM(stepnumber,1) = 28*S/pars.S0 ;
        workingstate.CO2ppm(stepnumber,1) = RCO2*280 ;
        workingstate.mrO2(stepnumber,1) = mrO2 ;
        workingstate.iceline(stepnumber,1) = iceline ;
        workingstate.T_gast(stepnumber,1) = TEMP_gast - 273 ;
        workingstate.ANOX(stepnumber,1) = ANOX ;
        workingstate.P(stepnumber,1) = P/pars.P0 ;
        workingstate.N(stepnumber,1) = N/pars.N0 ;
        workingstate.VEG(stepnumber,1) = VEG ;
        workingstate.time_myr(stepnumber,1) = t_geol ;
        workingstate.time(stepnumber,1) = t;
        workingstate.NPP(stepnumber,1) = sum( sum( NPP_past .* ( GRID_AREA_km2 .* 1e6 ) , 'omitnan') ) ; 
        workingstate.biomass(stepnumber,1) = sum( sum ( biomass_tot( end ), 'omitnan' ) ) ; 
        workingstate.habitablearea(stepnumber,1) = sum( sum( ~isnan( final_biomass_past ) .* ( GRID_AREA_km2 * 1e6 ) ) ) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Final actions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% output timestep if specified
if sensanal == 0
    if pars.telltime ==1
        if mod(stepnumber,pars.display_resolution) == 0 
            %%%% print model state to screen
            fprintf('Model step: %d \t', stepnumber); fprintf('time: %d \t', t_geol) ; fprintf('next keyframe: %d \n', next_stamp)
        end
    end
end

%%%% record current model step
stepnumber = stepnumber + 1 ;

%%%% option to bail out if model is running aground
% if stepnumber > pars.bailnumber
%    terminateExecution
% end

end

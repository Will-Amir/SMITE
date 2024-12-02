%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FLORA - Fast Land Occupancy Reaction Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BJW Mills 2021 (SCION)
%%%% b.mills@leeds.ac.uk
%%%% KG 2023 (FLORA)
%%%% k.gurung@leeds.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run = SCIFI_initialise_clean(runcontrol)

    %%%%%%% remove structures from pervious runs 
    clear stepnumber
    clear pars
    clear forcings
    clear workingstate
    clear switches
    clear state
    clear rawoutput
    clear options
    clear geoldata
    clear rawoutput
    clear resample
    clear biopars

    %%%%%%% set up global structures
    global stepnumber
    global pars
    global forcings
    global workingstate
    global state 
    global gridstate 
    global INTERPSTACK
    global sensanal
    global plotrun
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
    
    %%%%%%% global tuning variables
    global Gtune
    global Ctune
    global PYRtune
    global GYPtune
    global Atune
    global Otune
    global Stune
    
    %%%%%%% check for sensitivity analysis
    if runcontrol >= 1
        sensanal = 1 ;
        plotrun = 0 ;
        pars.telltime = 0 ;
    else
        sensanal = 0 ;
        plotrun = 1 ;
        pars.telltime = 1 ;
    end
    pars.runcontrol = runcontrol ; 
    
    %%%%%%% starting to load params
    if sensanal == 0 
        fprintf('setting parameters... \t')
        tic
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Flux values at present   %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% reductant input
    pars.k_reductant_input = 0.4e12 ; % schopf and klein 1992
    
    %%%%%%% org C cycle
    pars.k_locb = 2.5e12 ;
    pars.k_mocb = 2.5e12 ;
    pars.k_ocdeg = 1.25e12 ;

    %%%%%%% carb C cycle
    pars.k_ccdeg = 12e12 ;
    pars.k_carbw = 8e12 ;
    pars.k_sfw = 1.75e12 ;
    pars.k_mccb = pars.k_carbw + pars.k_ccdeg - pars.k_sfw ;
    pars.k_silw = pars.k_mccb - pars.k_carbw ;
    pars.basfrac = 0.3 ;
    pars.k_granw = pars.k_silw * (1-pars.basfrac) ;
    pars.k_basw = pars.k_silw * pars.basfrac ;

    %%%%%%% S cycle
    pars.k_mpsb = 0.7e12 ;
    pars.k_mgsb = 1e12 ;
    pars.k_pyrw = 7e11 ;
    pars.k_gypw = 1e12 ;
    pars.k_pyrdeg = 0 ; 
    pars.k_gypdeg = 0 ;
    
    %%%%%%% P cycle
    pars.k_capb = 2e10 ;
    pars.k_fepb = 1e10 ;
    pars.k_mopb = 1e10 ;
    pars.k_phosw = 4.25e10 ;
    pars.k_landfrac = 0.0588 ;
    %%%%%%% N cycle
    pars.k_nfix = 8.67e12 ;
    pars.k_denit = 4.3e12 ;

    %%%%%%% fluxes calculated for steady state
    pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg - pars.k_reductant_input ;

    %%%%%%% Sr cycle
    pars.k_Sr_sedw = 17e9 ;
    pars.k_Sr_mantle = 7.3e9 ;
    pars.k_Sr_silw = 13e9 ;
    pars.k_Sr_granw = pars.k_Sr_silw * (1 - pars.basfrac) ;
    pars.k_Sr_basw = pars.k_Sr_silw * pars.basfrac ;
    pars.total_Sr_removal = pars.k_Sr_granw + pars.k_Sr_basw + pars.k_Sr_sedw + pars.k_Sr_mantle ;
    pars.k_Sr_sfw = pars.total_Sr_removal * ( pars.k_sfw / (pars.k_sfw + pars.k_mccb) ) ;
    pars.k_Sr_sedb = pars.total_Sr_removal * ( pars.k_mccb / (pars.k_sfw + pars.k_mccb) ) ;
    pars.k_Sr_metam = 13e9 ;

    %%%%%%% others
    pars.k_oxfrac = 0.9975 ;
    Pconc0 = 2.2 ;
    Nconc0 = 30.9 ;
    pars.newp0 = 117 * min(Nconc0/16,Pconc0) ;
    % COPSE constant for calculating pO2 from normalised O2
    pars.copsek16 = 3.762 ;
    % oxidative weathering dependency on O2 concentration
    pars.a = 0.5 ;
    % marine organic carbon burial dependency on new production
    pars.b = 2 ; 
    % fire feedback
    pars.kfire= 3 ;

    % reservoir present day sizes (mol)
    pars.P0 = 3.1*10^15 ;
    pars.O0 = 3.7*10^19 ;
    pars.A0 = 3.193*10^18 ;
    pars.G0 = 1.25*10^21 ;
    pars.C0 = 5*10^21 ;
    pars.PYR0 = 1.8*10^20 ;
    pars.GYP0 = 2*10^20 ;
    pars.S0 = 4*10^19 ;
    pars.CAL0 = 1.397e19 ;
    pars.N0 = 4.35e16 ;
    pars.OSr0 = 1.2e17 ; % francois and walker 1992
    pars.SSr0 = 5e18 ;

    %%%%%% loading parameters for biomass calculations
    biopars.dt = 1 ;
    % Longitude and latitude 
    biopars.x_lon = 40 ; 
    biopars.y_lat = 48 ; 
    % Photosynthesis parameters (taken from LPJ model)
    biopars.t25 = 2600 ; 
    biopars.q_10t = 0.57 ; 
    biopars.v = 0.7 ; 
    biopars.kc25 = 30 ; 
    biopars.ko25 = 3e4 ; 
    biopars.q_10c = 2.1 ; 
    biopars.q_10o = 1.2 ; 
    biopars.s = (24/12) * 0.015 ; %0.08 in Github for C3 plants
    biopars.alpha = 0.08 ; 
    biopars.theta = 0.7 ; 
    biopars.lr_max = 0.75 ; 
    biopars.CN_leaf = 29 ;
    % Tissue respiration rate at 10 degree C
    biopars.r_tem = 0.055 * 365 ; % gC/gN/d -> gC/gN/year
    biopars.r_bor = 0.066 * 365 ;
    biopars.r_tro = 0.011 * 365 ; 
    % Growth respiraition
    biopars.R_growth = 0.25 ;     
    % Leaf longevity; ranges from 0.5 - 1 depending on type of plant
    biopars.life_leaf_tem = 0.75 ; 
    biopars.life_leaf_bor = 0.75 ;
    biopars.life_leaf_tro = 1 ; 
    % Minimum weathering
    biopars.minbiota = 0.32 ; 

    %%%%%%% finished loading params
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Forcings   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% starting to load forcings
    if sensanal == 0 
        fprintf('loading forcings... \t')
        tic
    end

    %%%%%%% load INTERPSTACK
    load( 'forcings/INTERPSTACK_oct_2021.mat' ) ;

    %%%%%%% relative contribution from latitude bands
    lat_areas = (cosd(INTERPSTACK.lat))' ;
    for n=1:48
        pars.rel_contrib(:,n) = lat_areas / mean(lat_areas) ;
    end
    
    %%%%%%% insolation parameter
    biopars.ins_present = 150 + 250 .* normpdf( INTERPSTACK.lat', 0 , 40 ) ./ normpdf( 0, 0, 40 ) .* ones( biopars.x_lon, biopars.y_lat ) ; 

    %%%%%%% load COPSE reloaded forcing set
    load( 'forcings/COPSE_forcings.mat' ) 
    %%% new BA 
    forcings.GR_BA = xlsread('forcings/GR_BA.xlsx','','','basic') ;
    forcings.GR_BA(:,1) = forcings.GR_BA(:,1)*1e6 ; %%% correct Myr
    %%% new GA
    forcings.newGA = xlsread('forcings/GA_revised.xlsx','','','basic') ;
    forcings.newGA(:,1) = forcings.newGA(:,1)*1e6 ; %%% correct Myr
    %%% degassing rate
    load('forcings/combined_D_force_oct_2021.mat') ;
    forcings.D_force_x = D_force_x ;
    forcings.D_force_mid = D_force_mid ;
    forcings.D_force_min = D_force_min ;
    forcings.D_force_max = D_force_max ;
    
    %%%%%%% load shoreline forcing
    load('forcings/shoreline.mat') ;
    forcings.shoreline_time = shoreline_time ;
    forcings.shoreline_relative = shoreline_relative ;

%     %%%%%%% load area
%     load('area.mat') ; 
    
    %%%%%%%% finished loading forcings
    if sensanal == 0 
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Sensitivity analysis   %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if sensanal == 1
        %%%% generate random number in [-1 +1]
        sensparams.randminusplus1 = 2*(0.5 - rand) ;
        sensparams.randminusplus2 = 2*(0.5 - rand) ;
        sensparams.randminusplus3 = 2*(0.5 - rand) ;
        sensparams.randminusplus4 = 2*(0.5 - rand) ;
        sensparams.randminusplus5 = 2*(0.5 - rand) ;
        sensparams.randminusplus6 = 2*(0.5 - rand) ;
        sensparams.randminusplus7 = 2*(0.5 - rand) ;    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%% run beginning
    if sensanal == 0 
        fprintf('Beginning run: \n')
    end

    %%%%%%% if no plot or sensitivity command set to single run
    if isempty(sensanal) == 1
        sensanal = 0 ;
    end
    if isempty(plotrun) == 1
        plotrun = 1 ;
    end

    %%%%%%% model timeframe in years (0 = present day)
    pars.whenstart = - 70e6 ;
    pars.whenend = -60e6 ;

    %%%%%%% impactor event properties (0 = present day)
    random_impactor_flag=0; % decides whether to randomly generate asteroids or if discerete asteroids should be used (0 is discrete, 1 is random)
    if random_impactor_flag==0
        asteroidtimes=[-67]; %Choose when in Myr to inject asteroids
        lats=[19]; %Choose asteroid latitudes
        longs=[13]; %Choose asteroid longitudes
        powers=[100]; %Choose asteroid powers
        SMITEflag=0;
        culledmaterial = ones(40,48) ;
        timetoinject=asteroidtimes(1);
    else
        numasteroids=1; %Choose number of random asteroids
        timetoinject=-60e6; %Choose when to inject asteroids
    end
    
    %%% Set seeding_primer to 0 for carrying forward biomass and 1 for reseeding
    seeding_primer=1;

    %%%%%%% setp up grid stamp times
    if runcontrol == -2
        pars.runstamps = 0 ;
    else
        pars.runstamps = INTERPSTACK.time( INTERPSTACK.time > ( pars.whenstart * 1e-6 ) ) ;
    end
    pars.next_gridstamp = pars.runstamps(1) ;
    pars.gridstamp_number = 1 ;
    pars.finishgrid = 0 ;
    
    %%%%%%% Show current timestep in command window? (1 = yes, 0 = no)
    pars.telltime = 1;

    %%%%%%% set number of model steps to take before beiling out
    pars.bailnumber = 1e5;

    %%%%%%% display every n model steps whilst running
    pars.display_resolution = 200 ;

    %%%%%%% set maximum step size for solver
    options = odeset('maxstep',1e6) ;

    %%%%%%% set stepnumber to 1
    stepnumber = 1 ;

    %%%%%%% set starting reservoir sizes 
    pars.pstart = pars.P0;
    pars.tempstart = 288;
    pars.CAL_start = pars.CAL0;
    pars.N_start = pars.N0;
    pars.OSr_start = pars.OSr0;
    pars.SSr_start = pars.SSr0;
    pars.delta_A_start = 0 ;
    pars.delta_S_start = 35 ;
    pars.delta_G_start = -27 ;
    pars.delta_C_start = -2 ;
    pars.delta_PYR_start = -5 ;
    pars.delta_GYP_start = 20 ;
    pars.delta_OSr_start = 0.708 ;
    pars.delta_SSr_start = 0.708 ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Initial parameter tuning  %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(Gtune) == 0
        pars.ostart = pars.O0 * abs( Otune )  ;
        pars.astart = pars.A0 * abs( Atune ) ;
        pars.sstart = pars.S0 * abs( Stune ) ;
        pars.gstart = pars.G0 * abs( Gtune ) ;
        pars.cstart = pars.C0 * abs( Ctune ) ;
        pars.pyrstart = pars.PYR0 * abs( PYRtune ) ;
        pars.gypstart = pars.GYP0 * abs( GYPtune ) ; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% if no tuning use previously tuned values
    if isempty(Gtune) == 1

        outputs = [ 0.55 0.96 1.2 1 0.1 0.05 3 ] ;  

        pars.gstart = pars.G0 * outputs(1) ;
        pars.cstart = pars.C0 * outputs(2) ;
        pars.pyrstart = pars.PYR0 * outputs(3) ;
        pars.gypstart = pars.GYP0 * outputs(4) ; 
        pars.ostart = pars.O0 * outputs(5)  ;
        pars.sstart = pars.S0 * outputs(6) ;
        pars.astart = pars.A0 * outputs(7) ;

    end

    %%%%%%% model start state
    pars.startstate(1) = pars.pstart ;
    pars.startstate(2) = pars.ostart ;
    pars.startstate(3) = pars.astart ;
    pars.startstate(4) = pars.sstart ;
    pars.startstate(5) = pars.gstart ;
    pars.startstate(6) = pars.cstart ;
    pars.startstate(7) = pars.pyrstart ;
    pars.startstate(8) = pars.gypstart ;
    pars.startstate(9) = pars.tempstart ;
    pars.startstate(10) = pars.CAL_start ;
    pars.startstate(11) = pars.N_start ;
    pars.startstate(12) = pars.gstart * pars.delta_G_start ;
    pars.startstate(13) = pars.cstart * pars.delta_C_start ;
    pars.startstate(14) = pars.pyrstart * pars.delta_PYR_start ;
    pars.startstate(15) = pars.gypstart * pars.delta_GYP_start ;
    pars.startstate(16) = pars.astart * pars.delta_A_start ;
    pars.startstate(17) = pars.sstart * pars.delta_S_start ;
    pars.startstate(18) = pars.OSr_start ;
    pars.startstate(19) = pars.OSr_start * pars.delta_OSr_start ;
    pars.startstate(20) = pars.SSr_start ;
    pars.startstate(21) = pars.SSr_start * pars.delta_SSr_start ;

    if seeding_primer == 0
        %%%%%%% initial FLORA running
        CO2ppm=18;
        timestart=18;
        mrO2 = ( pars.startstate(2)/pars.O0 )  /   ( (pars.startstate(2)/pars.O0)  + pars.copsek16 ) ;
        pO2 =mrO2*1000;
        timestartgeol=INTERPSTACK.time(timestart);
        Tair_past = INTERPSTACK.Tair(:,:,CO2ppm,timestart) ; 
        RUNOFF_past = INTERPSTACK.runoff(:,:,CO2ppm,timestart) ; 
        land_past = circshift(INTERPSTACK.land(:,:,timestart), [0 0]); 
        GRID_AREA_m2 = INTERPSTACK.gridarea(:,:)*1e6 ; 
        %setup
        dt = 1 ;
        % Longitude and latitude 
        x_lon = 48 ; 
        y_lat = 40 ;
        
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
        % Insolation parameter
        ins_present = 150 + 250 .* normpdf( INTERPSTACK.lat', 0 , 80 ) ./ normpdf( 0, 0, 80 ) .* ones( y_lat, x_lon) ; 
        
        
        %%%%%%% Insolation
        ins = ins_present - ( ins_present * 4.6/100 * ( abs( timestartgeol )/570 ) ) ; 
        
        %%%%%%% Ice (< -10 degC) or no runoff = no biomass 
        Tair_ice_past = Tair_past ; 
        Tair_ice_past( Tair_ice_past < -10 ) = NaN ;
        RUNOFF_past(RUNOFF_past == 0 | isnan(RUNOFF_past)) = NaN ; 
        RUNOFF_past = 1 - (1 ./ ( 1 + exp(0.005 .* (RUNOFF_past - 450)))) ;
        
        %%%%%%% Photosynthesis calculation for the past keyframe
        tf = t25 * ( q_10t .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 
        tstar = pO2 ./ ( 2 * tf ) ; 
        
        pi = v * CO2ppm ; 
        
        kc = kc25 * ( q_10c .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 
        ko = ko25 * ( q_10o .^ ( ( Tair_ice_past - 25 ) * 0.1 ) ) ; 
        
        c2 = ( pi - tstar ) ./ ( pi + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 
        
        ftemp_tem = normpdf( Tair_ice_past, 15,7 ) ; % 15, 7 ) ; % Temperate 
        ftemp_bor = normpdf( Tair_ice_past, 0, 20 ) ; %5, 10 ) ; % Boreal 
        ftemp_tro = normpdf( Tair_ice_past, 27, 7 ) ; % 27, 3 ) ;  % Tropical 
        
        c1_tem = alpha * ftemp_tem .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;
        c1_bor = alpha * ftemp_bor .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ; 
        c1_tro = alpha * ftemp_tro .* ( ( pi - tstar ) ./ ( pi + 2 * tstar ) ) ;
        
        sigma = ( 1 - ( ( c2 - s ) ./ ( c2 - theta * s ) ) ) .^ 0.5 ; 
        
        g_T_past = exp( 308.56 .* ( ( 1 / 56.02 ) - ( 1 ./ ( Tair_ice_past + 46.02 ) ) ) ) ;
        
        photosynth_tem_past = 10 .* 365 * ins .* ( c1_tem ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* RUNOFF_past ; 
        photosynth_bor_past = 10 .* 365 * ins .* ( c1_bor ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* RUNOFF_past ; 
        photosynth_tro_past = 10 .* 365 * ins .* ( c1_tro ./ c2 ) .* ( c2 - ( ( ( 2 * theta ) -1 ) * s ) - ( 2 .* ( c2 - ( theta .* s ) ) .* sigma ) ) .* RUNOFF_past ; 
        
        %%%%%%% Carbon in leaf allocation - first timestep
        C_leaf_tem_past = lr_max .* photosynth_tem_past ; 
        C_leaf_bor_past = lr_max .* photosynth_bor_past ;
        C_leaf_tro_past = lr_max .* photosynth_tro_past ; 
        
        %%%%%%% Biomass starting at homogenous values
        biomass_tem_past( :, :, 1 ) = 2.5e1 * land_past ; 
        biomass_tem_past( biomass_tem_past == 0 ) = NaN ; 
        biomass_bor_past( :, :, 1 ) = 2.5e1 * land_past ; 
        biomass_bor_past( biomass_bor_past == 0 ) = NaN ; 
        biomass_tro_past( :, :, 1 ) = 2.5e1 * land_past ; 
        biomass_tro_past( biomass_tro_past == 0 ) = NaN ; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%   Calculating Biomass   %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%% Stopping at <1% change in biomass
        
        biomass_change_final_step(1) = 999 ; 
        n = 1 ; 
        
        %%%%%%% Turnover changes between 8% and 20% depending on O2
        turnover = min( max( 0.0092 *( (pO2/1000) - 10 ), 0.08 ), 0.2 ) ; 
        
        while abs( biomass_change_final_step ) > 1 
        
            %%% Leaf respiration per plant type
            R_leaf_tem_past = r_tem * ( C_leaf_tem_past / CN_leaf ) .* g_T_past ;
            R_leaf_bor_past = r_bor * ( C_leaf_bor_past / CN_leaf ) .* g_T_past ;
            R_leaf_tro_past = r_tro * ( C_leaf_tro_past / CN_leaf ) .* g_T_past ;
            R_leaf_tem_past( R_leaf_tem_past <= 0 ) = 0 ; 
            R_leaf_bor_past( R_leaf_bor_past <= 0 ) = 0 ; 
            R_leaf_tro_past( R_leaf_tro_past <= 0 ) = 0 ; 
        
            %%% NPP
            NPP_tem_past = ( 1 - R_growth ) .* ( photosynth_tem_past - R_leaf_tem_past ) ; 
            NPP_bor_past = ( 1 - R_growth ) .* ( photosynth_bor_past - R_leaf_bor_past ) ; 
            NPP_tro_past = ( 1 - R_growth ) .* ( photosynth_tro_past - R_leaf_tro_past ) ; 
          
        
            %%% Biomass
            biomass_tem_past = biomass_tem_past + ( C_leaf_tem_past - turnover * biomass_tem_past) * dt ; 
            biomass_bor_past = biomass_bor_past + ( C_leaf_bor_past - turnover * biomass_bor_past ) * dt ; 
            biomass_tro_past = biomass_tro_past+ ( C_leaf_tro_past - turnover * biomass_tro_past ) * dt ; 
            
            final_biomass_past = max( biomass_tem_past, max( biomass_bor_past, biomass_tro_past ) ) ;
        
            %%% Carbon in leaf allocation (n+1)
            C_leaf_tem_past = ( C_leaf_tem_past .* ( 1 - life_leaf_tem ) ) + ( lr_max .* NPP_tem_past ) ;
            C_leaf_bor_past = ( C_leaf_bor_past .* ( 1 - life_leaf_bor ) ) + ( lr_max .* NPP_bor_past ) ;
            C_leaf_tro_past = ( C_leaf_tro_past .* ( 1 - life_leaf_tro ) ) + ( lr_max .* NPP_tro_past ) ;
        
        
            %%% Biomass total (gC/m2)
            biomass_tot( n + 1 ) = sum( sum( final_biomass_past .* ( GRID_AREA_m2), 'omitnan' ) ) ; 
            biomass_change_final_step = ( ( biomass_tot( n + 1 ) - biomass_tot( n ) ) / biomass_tot( n + 1 ) ) * 100 ; 
        
            n = n + 1 ; 
        end
        
        workingstate.biomass=final_biomass_past;
    end
    %%%%%%% note model start time
    tic

    %%%%%%% run the system 
    [rawoutput.T,rawoutput.Y] = ode15s(@SCIFI_equations2,[pars.whenstart pars.whenend],pars.startstate,options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% size of output 
    pars.output_length = length(rawoutput.T) ;
        
    if sensanal == 0
        %%%%%%% model finished output to screen
        fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
        toc
    end

    %%%%%%% print final model states using final state for each timepoint during integration
    
    if sensanal == 0
    fprintf('assembling state vectors... \t')
    %tic
    end

    %%%%%%% trecords is index of shared values between ode15s output T vector and model recorded workingstate t vector
    [sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

    %%%%%%% assemble output state vectors
    field_names = fieldnames(workingstate) ;
    for numfields = 1:length(field_names)
        eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
    end

    %%%%%%% save state
    run.state = state ;
    run.gridstate = gridstate ;
    run.pars = pars ;
    run.forcings = forcings ;
    
    if sensanal == 0
        %%% done message
        fprintf('Done: ')
        endtime = toc ;
        fprintf('time (s): %d \n', endtime )
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%   Plotting   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%% only plot if no tuning structure exists, only plot fluxes for quick runs
    if isempty(Gtune) == 1
        if plotrun == 1            
            if runcontrol>-1
                SCION_plot_worldgraphic
            end
            SCION_plot_fluxes_weathering
        end
    end

end

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
    clear options_long
    clear options_short
    clear geoldata
    clear rawoutput
    clear resample
    clear biopars
    clear biohistories
    clear biotimes

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
    global SMITE_asteroidparams
    global culledmaterial
    global seeding_primer
    global Crelease
    global growthmap

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
    
    %%%%%%% Atmospheric fraction of hydrosphere
    atfrac0 = 0.01614 ; % constant

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
    pars.atmoA0 = pars.A0*atfrac0 ;

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
    load( 'forcings/INTERPSTACK_CORE.mat' ) ;

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
    timetosimasteroid=5000; % sets how long after asteroid impact to simulate at high fidelity
    if random_impactor_flag==0
        SMITE_asteroidparams(1)=-67; %Choose when in Myr to inject asteroid
        SMITE_asteroidparams(2)=19; %Choose asteroid latitude
        SMITE_asteroidparams(3)=13; %Choose asteroid longitude
        SMITE_asteroidparams(4)=100; %Choose asteroid power
        SMITE_asteroidparams(5)=1e18; %Choose asteroid mass in grams
        culledmaterial = ones(40,48) ;
        growthmap=zeros(40,48);
    else
        numasteroids=1; %Choose number of random asteroids
        timetoinject=-60e6; %Choose when to inject asteroids
    end
    
    %%% Set seeding_primer to 0 for carrying forward biomass and 1 for reseeding
    seeding_primer=0;

    %%% Integration method, set to 0 for continuous integration and 1 for step
    stepflag=1;

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
    options_long = odeset('maxstep',1e6) ;
    options_short = odeset('MaxStep',1) ;


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
    pars.startstate(22) = pars.atmoA0 ;
    
    
    %%%%%%% note model start time
    tic

    %%%%%%% run the system 
    if stepflag==0
        [rawoutput.T,rawoutput.Y] = ode15s(@SCIFI_equations_standard,[pars.whenstart pars.whenend],pars.startstate,options_long);
    elseif stepflag==1
        timetostartinject=SMITE_asteroidparams(1)*1e6;
        timetoendinject=timetostartinject+timetosimasteroid;
        [rawoutput.T_init,rawoutput.Y_init] = ode15s(@SCIFI_equations_split, [pars.whenstart,timetostartinject],pars.startstate,options_long);
        SMITE
        %rawoutput.Y_init(end,22)=rawoutput.Y_init(end,22)+Crelease;
        [rawoutput.T_during,rawoutput.Y_during] = ode15s(@SCIFI_equations_regrowth, [timetostartinject,timetoendinject],rawoutput.Y_init(end,:),options_short);
        plot_script
        [rawoutput.T_post,rawoutput.Y_post] = ode15s(@SCIFI_equations_split, [timetoendinject,pars.whenend],rawoutput.Y_during(end,:),options_long);
        rawoutput.T=[rawoutput.T_init; rawoutput.T_during; rawoutput.T_post];
    end
    if runcontrol==-3
        FLORA_test
            run.singlerun = singlerun ; 
            run.pars = pars ;
            run.forcings = forcings ;
    end
    
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


accel = 1 ;
ratio_num = 3 ; 
max_thall = 600 ;
% light_fac = 0 ; 
% root_ratio = 0 ; 

if t_geol < -470 
%setting up environment for no plants    
    f_biota_past = biopars.minbiota .* RCO2 ^ 0.5 .* land_past  ;
    f_biota_future = biopars.minbiota .* RCO2 ^ 0.5 .* land_future ; 
    biomass_total = 0 ; 
    max_NPP_past = zeros(40,48) ; 
    final_biomass = zeros(40,48) ; 

    %%%%%%% Mass of biosphere
    VEG = 0 ; 
    turnover = 0 ; 
    ag_biomass = 0 .* land_past ; 
    bg_biomass = 0 .* land_past ; 
    biome = 0 .* land_past ; 
    maxIndex_past = zeros(40,48) ; 

 

else
    
    %%%%%%% Insolation
    ins = biopars.ins_present - ( biopars.ins_present * 4.6/100 * ( abs( t_geol )/570 ) ) ; 
    land_grid = land_past ; 
    land_grid(land_grid == 0 ) = NaN ; 
    land(:,:,1) = land_grid ; 
    land_grid = land_future ; 
    land_grid(land_grid == 0 ) = NaN ; 
    land(:,:,2) = land_grid ; 
    
    if t_geol >= -470 && t_geol < -400 
        %bryophytic roots 
         root_type = 1 ;
         bg_biomass = zeros(40,48,2,3,3) ; 
         ag_biomass = zeros(40,48,2,3,3) ;
         light_fac = 0.005 ; 
         root_ratio = 5 ;

     elseif t_geol >= -400 && t_geol < -350
         root_type = [1,2] ;
         bg_biomass = zeros(40,48,2,3,3) ;
         ag_biomass = zeros(40,48,2,3,3) ; 
         light_fac = 0.005 ; 
         root_ratio = 5 ; 

     elseif t_geol > -350  
         root_type = [1,2,3] ; 
         light_fac = 0.005 ; 
         root_ratio = 5 ; 
    end
   % 3D = time, 4D = biome, 5D = root type
   % 1 = past, 2 = future
   % 1 = boreal, 2 = temperate, 3 = tropical
   % 1 = rhizoid, 2 = tracheophytes, 3 = seed plants

    %%%%%%% Ice (< -10 degC) or no runoff = no biomass 
    Tair_ice_past = Tair_past .* land(:,:,1); 
    Tair_ice_future = Tair_future .* land(:,:,2) ; 
    Tair_ice_past( TOPO_past > 4000 ) = NaN ; 
    Tair_ice_future( TOPO_future > 4000 ) = NaN ;
    Tair_ice_past( Tair_ice_past < -20 ) = NaN ; 
    Tair_ice_future( Tair_ice_future < -20 ) = NaN ; 
    Tair_ice(:,:,1,:,root_type) = repmat(Tair_ice_past,1,1,1,3,length(root_type)) ; 
    Tair_ice(:,:,2,:,root_type) = repmat(Tair_ice_future,1,1,1,3,length(root_type)) ; 

    water_depth(:,:,1,:,root_type) = repmat( 83 .* exp( -0.01 .* RUNOFF_past ), 1, 1, 1, 3, length(root_type) ) ; %0.01
    water_depth(:,:,2,:,root_type) = repmat( 83 .* exp( -0.01 .* RUNOFF_future ), 1, 1, 1, 3, length(root_type) ) ; 
    %Rooting depth from Canadell et al. (1996)
    rel_rooting_depth(:,:,:,1,root_type) = biopars.rooting_depth(:,:,:,1,root_type) / 4 ; %boreal
    rel_rooting_depth(:,:,:,2,root_type) = biopars.rooting_depth(:,:,:,2,root_type) / 6 ; %temperate
    rel_rooting_depth(:,:,:,3,root_type) = biopars.rooting_depth(:,:,:,3,root_type) / 8 ; %tropical
    
    water_avail(:,:,1,:,root_type) = 1 - (1 ./ ( 1 + exp(1 .* (biopars.rooting_depth(:,:,1,:,root_type) - water_depth(:,:,1,:,:)))) .* land(:,:,1,:,:)) ;  
    water_avail(:,:,2,:,root_type) = 1 - (1 ./ ( 1 + exp(1 .* (biopars.rooting_depth(:,:,2,:,root_type) - water_depth(:,:,2,:,:)))) .* land(:,:,2,:,:)) ; 
    water_avail(water_avail == 0) = NaN ; 

    %%%%%%% Photosynthesis calculation for the past and future keyframe
    tf = biopars.t25 * ( biopars.q_10t .^ ( ( Tair_ice - 25 ) .* 0.1 ) ) ;
    tstar = pO2 ./ ( 2 .* tf ) ; 
    
    intra_pp = biopars.v * CO2ppm ; 
   
    if t_geol < -400
        kc = biopars.kc25_lp .* ( biopars.q_10c_lp .^ ( ( Tair_ice(:,:,:,:,1) - 25 ) .* 0.1 ) ) ;
        ko = biopars.ko25_lp .* ( biopars.q_10o_lp .^ ( ( Tair_ice( :, :, :, :, 1 ) - 25 ) .* 0.1 ) ) ; 
    else
        root_type(1) = [] ; 
        kc_lp = biopars.kc25_lp .* ( biopars.q_10c_lp .^ ( ( Tair_ice(:,:,:,:,1) - 25 ) .* 0.1 ) ) ;
        kc_hp = biopars.kc25 .* ( biopars.q_10c .^ ( ( Tair_ice(:,:,:,:,root_type) - 25 ) .* 0.1 ) ) ; 
        kc = cat(5, kc_lp, kc_hp) ; 
        ko_lp = biopars.ko25_lp .* ( biopars.q_10o_lp .^ ( ( Tair_ice( :, :, :, :, 1 ) - 25 ) .* 0.1 ) ) ; 
        ko_hp = biopars.ko25 .* ( biopars.q_10o .^ ( ( Tair_ice( :, :, :, :, root_type) - 25 ) .* 0.1 ) ) ; 
        ko = cat(5, ko_lp, ko_hp) ; 
    end
    c2 = ( intra_pp - tstar ) ./ ( intra_pp + kc .* ( 1 + ( pO2 ./ ko ) ) ) ; 
    
    clear ftemp
    ftemp(:,:,:,1,:) = normpdf( Tair_ice(:,:,:,1,:), 1, 10 ) ; %5, 10 ) ; % Boreal  
    ftemp(:,:,:,2,:) = normpdf( Tair_ice(:,:,:,2,:), 10, 7 ) ; % 15, 7 ) ; % Temperate 
    ftemp(:,:,:,3,:) = normpdf( Tair_ice(:,:,:,3,:), 27, 3 ) ; % 27, 3 ) ;  % Tropical

    c1 = biopars.alpha .* ftemp .* ( ( intra_pp - tstar ) ./ ( intra_pp + 2 .* tstar ) ) ;

    sigma = min(0, ( 1 - ( ( c2 - biopars.s ) ./ ( c2 - biopars.theta .* biopars.s ) ) ) .^ 0.5 ) ; 
    Arr_eq = exp( 308.56 .* ( ( 1 / 56.02 ) - ( 1 ./ ( Tair_ice + 46.02 ) ) ) ) ;
    %Soil estimated temperature
    T_soil = Tair_ice + exp( -biopars.z / biopars.d ) .* sin( biopars.omega - ( biopars.z / biopars.d ) ) ;
    gT_soil = exp( 308.56 .*( ( 1 / 56.02 ) - ( 1 ./ ( T_soil + 46.02 ) ) ) ) ;
    
    %metabolic activity limiation from Porada
    met_act = min(1, water_avail(:,:,:,:,1) ./ 0.3 ) ;  

    LAI = ones(40,48,2,3,3) ;
    L = exp(-light_fac .* LAI) ; 

        %%%%%%% Carbon allocation - first timestep
    Croot_ratio = rel_rooting_depth(:,:,:,:,root_type) .* root_ratio .* 0.3 .* ( L(:,:,:,:,root_type) ./ (L(:,:,:,:,root_type)  + 2 .* water_avail(:,:,:,:,root_type)) ) ;
    Cstem_ratio = ratio_num .* 0.3 .* ( water_avail(:,:,:,:,root_type) ./ (2 .* L(:,:,:,:,root_type)  + water_avail(:,:,:,:,root_type) ) ) ; 
    Cleaf_ratio = 1 - ( Croot_ratio + Cstem_ratio ) ;

    if t_geol<-400
        photosynth =  met_act .* water_avail(:,:,:,:,1) .* land .* 365 .* ins .* ( c1(:,:,:,:,1) ./ c2(:,:,:,:,1) ) .* ( c2(:,:,:,:,1) - ( ( ( 2 .* biopars.theta ) -1 ) .* biopars.s ) - ( 2 .* ( c2(:,:,:,:,1) - ( biopars.theta .* biopars.s ) ) .* sigma(:,:,:,:,1) ) ) ; 
        
        C_thalloid = min(max_thall, photosynth,"includenan") .* 1 ; %LAI of 1  
        clear ag_biomass bg_biomass
        ag_biomass( :, :,:,:, : ) = zeros(40,48, 2,3,1) .* land ;  
        bg_biomass( :, :,:,:, : ) = zeros(40,48, 2,3,1) .* land ; 
    else
        photosynth_lp =  met_act .* water_avail(:,:,:,:,1) .* land .* 365 .* ins .* ( c1(:,:,:,:,1) ./ c2(:,:,:,:,1) ) .* ( c2(:,:,:,:,1) - ( ( ( 2 .* biopars.theta ) -1 ) .* biopars.s ) - ( 2 .* ( c2(:,:,:,:,1) - ( biopars.theta .* biopars.s ) ) .* sigma(:,:,:,:,1) ) ) ; 
        photosynth_hp =  water_avail(:,:,:,:,root_type) .* land .* 365 .* ins .* ( c1(:,:,:,:,root_type) ./ c2(:,:,:,:,root_type) ) .* ( c2(:,:,:,:,root_type) - ( ( ( 2 .* biopars.theta ) -1 ) .* biopars.s ) - ( 2 .* ( c2(:,:,:,:,root_type) - ( biopars.theta .* biopars.s ) ) .* sigma(:,:,:,:,root_type) ) ) ; 
        photosynth = cat(5,photosynth_lp, photosynth_hp) ; 

        %Leaf
        C_leaf = Cleaf_ratio .* photosynth(:,:,:,:,root_type) .* LAI(:,:,:,:,root_type) ; 
     
        %Root
%         C_root = rel_rooting_depth(:,:,:,:,root_type) .* Croot_ratio .* photosynth(:,:,:,:,root_type) .* LAI(:,:,:,:,root_type); 
        C_root = Croot_ratio .* photosynth(:,:,:,:,root_type) .* LAI(:,:,:,:,root_type); 
        %Stem
        C_stem = Cstem_ratio .* photosynth(:,:,:,:,root_type) .* LAI(:,:,:,:,root_type); 
    
        %Thalloid only
        C_thalloid = min(max_thall, photosynth(:,:,:,:,1),"includenan") .* LAI(:,:,:,:,1) ;
        
        clear ag_biomass bg_biomass
        ag_biomass( :, :,:,:, : ) = zeros(40,48, 2,3,length(root_type)+1) .* land ;  
        bg_biomass( :, :,:,:, : ) = zeros(40,48, 2,3, length(root_type)+1) .* land ; 

        
    end



    
    
    %%%%%%% Biomass starting at homogenous values

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%   Calculating Biomass   %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% Stopping at <1% change in biomass
    biomass_change_final_step(1) = 999 ; 
    n = 1 ; 
    
    %%%%%%% Turnover changes between 8% and 20% depending on O2
 %   turnover = min( max( 0.0092 *( mrO2 * 100 - 10 ), 0.08 ), 1 ) ; 
    turnover = sigmf(mrO2, [4.5 2]) * 1000 + 0.03 ; 
%     turnover = min( max( 0.03 * (mrO2 * 100 - 10 ), 0.08 ), 1 ) ; 
    [rows,cols] = meshgrid(1:48, 1:40) ; 
    
    if runcontrol == -3
        perc_change = 0.01 ; 
    else
        perc_change = 1; 
    end

    while abs( biomass_change_final_step ) > perc_change   

                if t_geol < -400
                    total_resp = biopars.r(:,:,:,:,1).^( (Tair_ice(:,:,:,:,1) - 9.85) ./ (9.85 + 0.15 ) ) .* C_thalloid .* met_act ; 
                    NPP = ( 1 - biopars.R_growth ) .* ( photosynth - total_resp ) ; 
                    ag_biomass(:,:,:,:,1) = ag_biomass(:,:,:,:,1) + ( ( C_thalloid ) - turnover .* ag_biomass(:,:,:,:,1) ) .* biopars.dt ; 
                    bg_biomass(:,:,:,:,1) = bg_biomass(:,:,:,:,1) - turnover .* bg_biomass(:,:,:,:,1) .* biopars.dt ;
                    temporary_grid = NPP > 0 ; 
            C_thalloid = min(C_thalloid .* (1 - biopars.life_thall) + (NPP(:,:,:,:,1) .* temporary_grid(:,:,:,:,1)) .* LAI(:,:,:,:,1) ,max_thall,'includenan') ; 

        else    
     %%% Leaf respiration per plant type
        %Temperate
            R_leaf = biopars.r(:,:,:,:,root_type) .* ( C_leaf ./ biopars.CN_leaf ) .* Arr_eq(:,:,:,:,root_type) ;
            R_root = biopars.r(:,:,:,:,root_type) .* ( C_root ./ biopars.CN_leaf ) .* gT_soil(:,:,:,:,root_type) ; 
            R_stem = biopars.r(:,:,:,:,root_type) .* ( C_stem ./ biopars.CN_sapwood ) .* Arr_eq(:,:,:,:,root_type) ; 
            R_thall = biopars.r(:,:,:,:,1).^( (Tair_ice(:,:,:,:,1) - 9.85) ./ (9.85 + 0.15 ) ) .* C_thalloid .* met_act ; 
            total_resp = cat(5, R_thall, R_leaf + R_root + R_stem) ; 
            %%% Update LAI parameters
            %Temperate
            %LAI_bor = sigmf(LAI_bor,[0.01 750]) .* 2 + 3 ;  
            LAI(:,:,:,1,root_type) = 1 ./ ( 1 + exp( -0.01 .* ( (C_leaf(:,:,:,1,:) + C_stem(:,:,:,1,:)) - 750 ) ) ) .* 2 + 3 ; 
            %LAI_tem = sigmf(LAI_tem,[0.01 750]) .* 2 + 2 ;
            LAI(:,:,:,2,root_type) = 1 ./ ( 1 + exp( -0.01 .* ( (C_leaf(:,:,:,2,:) + C_stem(:,:,:,2,:)) - 750 ) ) ) .* 2 + 2 ; 
            %LAI_tro = sigmf(LAI_tro,[0.005 1500]) .* 2 + 4 ; 
            LAI(:,:,:,3,root_type) = 1 ./ ( 1 + exp( -0.01 .* ( (C_leaf(:,:,:,3,:) + C_stem(:,:,:,3,:)) - 1500 ) ) ) .* 2 + 4 ; 
%             LAI(:,:,:,:,1) = 1./ (1 + exp( - 0.01 .* ((C_thalloid - 750 ) ) ) ) .* 2 ; 
            L = exp( -light_fac * LAI) ;
            %%% NPP
            NPP = ( 1 - biopars.R_growth ) .* ( photosynth - total_resp ) ; 
        
       
            %%% Biomass
            ag_biomass(:,:,:,:,1) = ag_biomass(:,:,:,:,1) + ( C_thalloid - turnover .* ag_biomass(:,:,:,:,1) ) .* biopars.dt ; 
            ag_biomass(:,:,:,:,root_type) = ag_biomass(:,:,:,:,root_type) + ( ( C_leaf + C_stem ) - turnover .* ag_biomass(:,:,:,:,root_type) ) .* biopars.dt ;  
            bg_biomass(:,:,:,:,1) = bg_biomass(:,:,:,:,1) - turnover .* bg_biomass(:,:,:,:,1) .* biopars.dt ; 
            bg_biomass(:,:,:,:,root_type) = bg_biomass(:,:,:,:,root_type) + ( C_root - turnover.* bg_biomass(:,:,:,:,root_type) ) .* biopars.dt ; 
        
   
        %%% C allocation %%%
            temporary_grid = NPP > 0 ; 
            C_leaf = C_leaf .* ( 1 - biopars.life_leaf(:,:,:,:,length(root_type)) ) + ( NPP(:,:,:,:,root_type) .* temporary_grid(:,:,:,:,root_type) .* Cleaf_ratio .* (LAI(:,:,:,:,root_type) ./ biopars.maxLAI(:,:,:,:,length(root_type))) ) ; 
%             C_root = C_root .* ( 1 - biopars.life_root(:,:,:,:,length(root_type)) ) + ( NPP(:,:,:,:,root_type) .* temporary_grid(:,:,:,:,root_type) .* Croot_ratio .* rel_rooting_depth(:,:,:,:,root_type) .* (LAI(:,:,:,:,root_type) ./ biopars.maxLAI(:,:,:,:,length(root_type))) )  ; 
            C_root = C_root .* ( 1 - biopars.life_root(:,:,:,:,length(root_type)) ) + ( NPP(:,:,:,:,root_type) .* temporary_grid(:,:,:,:,root_type) .* Croot_ratio .* (LAI(:,:,:,:,root_type) ./ biopars.maxLAI(:,:,:,:,length(root_type))) )  ; 
            C_stem = C_stem .* ( 1 - biopars.life_sapwood(:,:,:,:,length(root_type)) ) + ( NPP(:,:,:,:,root_type) .* temporary_grid(:,:,:,:,root_type) .* Cstem_ratio .* (LAI(:,:,:,:,root_type) ./ biopars.maxLAI(:,:,:,:,length(root_type)))  ) ; 
            C_thalloid = min(C_thalloid .* (1 - biopars.life_thall) + (NPP(:,:,:,:,1) .* temporary_grid(:,:,:,:,1)) .* LAI(:,:,:,:,1) ,max_thall,'includenan') ; 
            %Update LAI and NPP allocation ratios for next step
            %Temperate
            Croot_ratio = rel_rooting_depth(:,:,:,:,root_type) .* root_ratio .* 0.3 .* ( L(:,:,:,:,root_type) ./ (L(:,:,:,:,root_type) + 2 .* water_avail(:,:,:,:,root_type) ) ) ; 
            Cstem_ratio = ratio_num .* 0.3 .* ( water_avail(:,:,:,:,root_type) ./ (2 .* L(:,:,:,:,root_type) + water_avail(:,:,:,:,root_type) ) ) ; 
            Cleaf_ratio = 1 - ( Croot_ratio + Cstem_ratio ) ; 
        end

        %%%% reshape T and roots dependencies into single stack
        ag_biomass_past = cat( 3, NaN(40,48), reshape( ag_biomass(:,:,1,:,:), 40,48,[] ) ) ; %reshaping puts it in order; eg: [1 2 3; 4 5 6] becomes [1 2 3 4 5 6]; %first sheet all NaNs 
        ag_biomass_future = cat( 3, NaN(40,48), reshape( ag_biomass(:,:,2,:,:), 40,48,[] ) ) ;
        [max_ag_past, maxIndex_past] = max(ag_biomass_past,[],3) ; %takes the max out of all sheets in 3rd dimension
        [max_ag_future, maxIndex_future] = max(ag_biomass_future,[],3) ; 

        %%%% find corresponding bg biomass
        bg_biomass_past = cat( 3, NaN(40,48), reshape( bg_biomass(:,:,1,:,:), 40,48,[] ) ) ;
        bg_biomass_future = cat( 3, NaN(40,48), reshape( bg_biomass(:,:,2,:,:), 40,48,[] ) ) ;
%         bg_biomass_past_reshape = reshape(bg_biomass_past, 40,48,[]) ;
%         bg_biomass_future_reshape = reshape(bg_biomass_future, 40,48,[]) ;
        
        linearIndices = sub2ind(size(bg_biomass_past),  cols(:), rows(:), maxIndex_past(:)) ; 
        max_bg_past = reshape(bg_biomass_past(linearIndices), [40, 48, 1]) ; 
        linearIndices = sub2ind(size(bg_biomass_past),  cols(:), rows(:), maxIndex_future(:)) ; 
        max_bg_future = reshape(bg_biomass_future(linearIndices), [40,48,1]) ; 
  
        final_biomass(:,:,1) = max_ag_past + max_bg_past ; 
        final_biomass(:,:,2) = max_ag_future + max_bg_future ; 
        %%% Biomass total (gC/m2)

        biomass_total( n + 1 ) = sum(sum(final_biomass(:,:,1) .* GRID_AREA_km2 * 1e6, 'omitnan')) * contribution_past + sum(sum(final_biomass(:,:,2) .* GRID_AREA_km2 * 1e6,'omitnan')) * contribution_future ; 
        biomass_change_final_step = ( ( biomass_total( n + 1 ) - biomass_total( n ) ) / biomass_total( n + 1 ) ) * 100 ; 
  
        if runcontrol == -3 
            if t_geol < -400 
                singlerun.biomass_total(n) = biomass_total(n) ;
                singlerun.final_biomass(:,:,n) = final_biomass(:,:,1) ; 
                singlerun.biome(:,:,n) = maxIndex_past ; 
                singlerun.ag_biomass(:,:,n) = max_ag_past ; 
                singlerun.bg_biomass(:,:,n) = max_bg_past ;
                singlerun.thall(:,:,n) = C_thalloid(:,:,1,3,1) ;
            else
                singlerun.biomass_total(n) = biomass_total(n) ;
                singlerun.final_biomass(:,:,n) = final_biomass(:,:,1) ; 
                singlerun.biome(:,:,n) = maxIndex_past ; 
                singlerun.ag_biomass(:,:,n) = max_ag_past ; 
                singlerun.bg_biomass(:,:,n) = max_bg_past ; 
                singlerun.root(:,:,n) = C_root(:,:,1,3,2) ;
                singlerun.stem(:,:,n) = C_stem(:,:,1,3,2) ; 
                singlerun.leaf(:,:,n) = C_leaf(:,:,1,3,2) ; 
            end 
        end



        n = n + 1 ; 
%         fprintf('n= %d \n' , n )
        accel = accel + 10;


    end
 

 if runcontrol == -3
     n
     if t_geol >-400
     figure
     subplot(4,5,1)
     imagesc(Cstem_ratio(:,:,1,3,2))
     colorbar
     title('Stem ratio (Tro;Full)')
     subplot(4,5,2)
     imagesc(Cleaf_ratio(:,:,1,3,2))
     colorbar
     title('Leaf ratio')
     subplot(4,5,3)
     imagesc(Croot_ratio(:,:,1,3,2))
     colorbar
     title('Root ratio')
     subplot(4,5,4)
     imagesc(C_stem(:,:,1,3,2))
     colorbar
     title('C stem')
     subplot(4,5,5)
     imagesc(C_leaf(:,:,1,3,2))
     colorbar
     title('C leaf')
     subplot(4,5,6)
     imagesc(C_root(:,:,1,3,2))
     colorbar
     title('C root')
     end
 end
%%%% Finding NPP for the used root/temp combination

    NPP_past = cat(3, NaN(40,48), reshape(NPP(:,:,1,:,:), 40,48,[] ) ) ; 
    NPP_future = cat(3, NaN(40,48), reshape( NPP(:,:,2,:,:), 40,48,[] ) ) ; 
    max_NPP_future = reshape(NPP_future(linearIndices), [40,48,1]) ; 
    linearIndices = sub2ind(size(NPP_past), cols(:), rows(:), maxIndex_past(:)) ;
    max_NPP_past = reshape(NPP_past(linearIndices), [40,48,1]) ; 
%     C_leaf_past = reshape(C_thalloid(:,:,1,:),C_leaf(:,:,1,:,:), 40,48,[]) ;
%     linearIndicies = sub2ind(size(C_leaf_past), cols(:), rows(:), maxIndex_past(:)) ; 
%     max_Cleaf = reshape(C_leaf_past(linearIndices), [40,48,1]) ;
%     C_stem_past = reshpae(C_stem(:,:,1,:,:), 40,48,[]) ;
%     linearIndicies = sub2ind(size(C_stem_past), cols(:), rows(:), maxIndex_past(:)) ; 
%     max_Cstem = reshape(C_stem_past(linearIndices), [40,48,1]) ;
%     C_root_past = reshpae(C_root(:,:,1,:,:), 40,48,[]) ;
%     linearIndicies = sub2ind(size(C_root_past), cols(:), rows(:), maxIndex_past(:)) ; 
%     max_Croot = reshape(C_root_past(linearIndices), [40,48,1]) ;

    %%%%%%% Calculating effect of biomass on weathering
    
% 
%     total_biomass_tem_past(total_biomass_tem_past<0.1) = NaN ; 
%     alloc_grid = final_biomass_past == total_biomass_tem_past ;
%     ag_biomass = alloc_grid .* ag_biomass_tem_past(:,:,end) ; 
%     bg_biomass = alloc_grid .* bg_biomass_tem_past(:,:,end) ; 
%     NPP_past = alloc_grid .* NPP_tem(:,:,1) ; 
%     clear LAI
%     LAI(:,:,1) = alloc_grid .* LAI_tem(:,:,1) ; 
%     biome = alloc_grid .* 1 ; 
%     total_biomass_tro_past(total_biomass_tro_past<0.1) = NaN ; 
%     alloc_grid = final_biomass_past == total_biomass_tro_past ;
%     ag_biomass = ag_biomass + (alloc_grid .* ag_biomass_tro_past(:,:,end)) ; 
%     bg_biomass = bg_biomass + (alloc_grid .* bg_biomass_tro_past(:,:,end)) ; 
%     NPP_past = NPP_past + (alloc_grid .* NPP_tro(:,:,1)) ; 
%     LAI(:,:,1) = LAI(:,:,1) + (alloc_grid .* LAI_tro(:,:,1)) ;
%     biome = biome + alloc_grid .* 3 ; 
%     total_biomass_bor_past(total_biomass_bor_past<0.1) = NaN ; 
%     alloc_grid = final_biomass_past == total_biomass_bor_past ;
%     ag_biomass = ag_biomass + (alloc_grid .* ag_biomass_bor_past(:,:,end)) ; 
%     bg_biomass = bg_biomass + (alloc_grid .* bg_biomass_bor_past(:,:,end)) ; 
%     NPP_past = NPP_past + (alloc_grid .* NPP_bor(:,:,1)) ; 
%     LAI(:,:,1) = LAI(:,:,1) + (alloc_grid .* LAI_bor(:,:,1)) ;
%     biome = biome + alloc_grid .* 2 ; 
% 
%     total_biomass_tem_future(total_biomass_tem_future<0.1) = NaN ; 
%     alloc_grid = final_biomass_future == total_biomass_tem_future ; 
%     NPP_future = alloc_grid .* NPP_tem(:,:,2) ; 
%     LAI(:,:,2) = alloc_grid .* LAI_tem(:,:,2) ; 
%     total_biomass_tro_future(total_biomass_tro_future<0.1) = NaN ; 
%     alloc_grid = final_biomass_future == total_biomass_tro_future ; 
%     NPP_future = NPP_future + (alloc_grid .* NPP_tro(:,:,2)) ; 
%     LAI(:,:,2) = LAI(:,:,2) + (alloc_grid .* LAI_tro(:,:,2)) ; 
%     total_biomass_bor_future(total_biomass_bor_future<0.1) = NaN ; 
%     alloc_grid = final_biomass_future == total_biomass_bor_future ; 
%     NPP_future = NPP_future + (alloc_grid .* NPP_bor(:,:,2)) ;
%     LAI(:,:,2) = LAI(:,:,2) + (alloc_grid .* LAI_bor(:,:,2)) ; 

    f_biota_past = 0.001 * max_NPP_past + ( biopars.minbiota * RCO2 ^ 0.6 ) ; 
    f_biota_future = 0.001 * max_NPP_future + ( biopars.minbiota * RCO2 ^ 0.6 ) ; 
    
    f_biota_past( land_past == 1 & isnan( f_biota_past ) ) = biopars.minbiota .* RCO2 ^ 0.6 ;
    f_biota_future( land_future == 1 & isnan( f_biota_future ) ) = biopars.minbiota .* RCO2 ^ 0.6 ; 
    
    %%%%%%% Mass of biosphere
    VEG = biomass_total( end ) / 4.5e17 ; %2e18 ; 
    %best one with 3.5e18, 0.6OC, 0.92C

%     if runcontrol == -3
%         singlerun.Cleaf = max_Cleaf ;
%         singlerun.Croot = max_Croot ; 
%         singlerun.Cstem = max_Cstem ; 
% %         singlerun.biome = 
% %         singlerun.
%     end
end

% a = [1 2 3 ; 13 23 18; 19 27 12] ; 
% b = [20 24 11; 4 5 6; 14 25 17];
% c = [10 26 16; 21 22 15; 7 8 9] ; 
% d = [100 101 102; 125 110 111; 116 120 109] ; 
% e = [121 117 112 ; 103 108 119; 124 115 104] ; 
% f = [113 106 114; 105 122 107; 123 118 126] ; 
% g = cat(3, a, b, c) ;
% h = cat(3, d, e,f ) ;
% i = cat(4, g, h) ;
% test_reshape = reshape(i, 3,3,[]) ;
% j = [1 2 3 ; 4 5 6 ; 7 8 9] ; 
% k = [10 11 12 ; 13 14 15 ; 16 17 18] ; 
% l = [19 20 21 ; 22 23 24 ; 25 26 27] ; 
% m = [28 29 30 ; 31 32 33 ; 34 35 36] ; 
% n = [37 38 39 ; 40 41 42 ; 43 44 45] ; 
% o = [46 47 48 ; 49 50 51 ; 52 53 54] ; 
% p = cat(3, j,k,l) ; 
% q = cat(3, m, n, o) ; 
% r = cat(4, p,q) ; 
% test_reshape_2 = reshape(r, 3,3,[]) ;
% [rows, cols] = meshgrid(1:3, 1:3) ; 
% linearIndices = sub2ind(size(test_reshape_2),  cols(:), rows(:), testindex(:)) ; 
% resultArray = test_reshape_2(linearIndices) ; 
% resultArray = reshape(resultArray, [3 3 1]) ;

% for i = 1 : 6
%     im_land = gridstate.land(:,:,i+2) ; 
%     im_land(im_land == 0) = NaN ; 
%     im_b = gridstate.BIOMASS_tot_grid(:,:,i+2) ; 
%     im_b(isnan(im_b) == 1 ) = 1 ; 
%     subplot(3,2,i)
%     h = pcolor(flipud(im_land .* im_b)) ; 
%     set(h, 'EdgeColor', 'none')
%     title(gridstate.time_myr(i+2))
%     colorbar
%     caxis([0 5e4])
% 
% end
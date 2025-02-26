%%%%%%%%%%%%%%%%SMITE
clear toburn
clear biomasstocarry
global biopars
global culledmaterial
global INTERPSTACK
global toburn
global Crelease
global biomasstocarry
global growthmap
global SMITE_asteroidparams
global atmoCO2

GRID_AREA_km2=INTERPSTACK.gridarea;
key_time = max ( INTERPSTACK.time( (INTERPSTACK.time - SMITE_asteroidparams(1)) <= 0 ) ) ;
key_index = find( INTERPSTACK.time == key_time ) ;
key_CO2 = max (INTERPSTACK.CO2( ( INTERPSTACK.CO2-atmoCO2)<=0));
CO2_index = find(INTERPSTACK.CO2==key_CO2);
plantgridarea=GRID_AREA_km2;
land_past = INTERPSTACK.land(:,:,key_index) ; 
potentialfires = INTERPSTACK.wildfires(:,:) ;
topologymap=INTERPSTACK.topo(:,:,key_index);
runoffmap=INTERPSTACK.runoff(:,:,CO2_index,key_index);
struck=INTERPSTACK.struck;
if random_impactor_flag==1
    lats = zeros(1,length(asteroidtimes)) ;
    longs = zeros(1,length(asteroidtimes)) ;
    powers = zeros(1,length(asteroidtimes)) ;
    %Generates asteroids
    for i = 1:length(asteroidtimes)
        lats(i) = randi(biopars.y_lat) ;
        longs(i) = randi(biopars.x_lon) ;
        powers(i) = randi(100) ;
        struck(lats(i),longs(i))=1;
    end
else
    asteroidtime=SMITE_asteroidparams(1);
    lat=SMITE_asteroidparams(2);
    long=SMITE_asteroidparams(3);
    power=SMITE_asteroidparams(4);
    mass=SMITE_asteroidparams(5);
    struck(lat,long)=1;
end
struck=circshift(struck, [0 29]);
%Generates surrounding cells for wildires
struckcells = find(struck==1) ;
struckcellsx = floor(struckcells/biopars.x_lon) ;
struckcellsy = rem(struckcells,biopars.x_lon) ;
for j = 1:length(asteroidtime)
    if struckcellsy(j)==biopars.y_lat
        toaddy=[39,40];
    elseif struckcellsy(j)==0
        toaddy=[1,2];
    elseif struckcellsy(j)==1
        toaddy=[1,2,3];
    else
        toaddy=[struckcellsy(j)-3,struckcellsy(j)-2,struckcellsy(j)-1,struckcellsy(j),struckcellsy(j)+1,struckcellsy(j)+2,struckcellsy(j)+3];
    end
    if struckcellsx(j)==29
        toaddx=[28,29,27,26];
    elseif struckcellsx(j)==28
        toaddx=[26,27,28,29];
    else
        toaddx=[struckcellsx(j)-3,struckcellsx(j)-2,struckcellsx(j)-1,struckcellsx(j),struckcellsx(j)+1,struckcellsx(j)+2,struckcellsx(j)+3];
    end
    if power>=98
        potentialfires(:,:) = 1;
    else
        potentialfires(toaddy,toaddx) = 1;
    end
end
%Checks that wildfires are on land
truefires=potentialfires ;
truefires(potentialfires==1 & land_past==1)=1 ;
truefires(potentialfires==0 & land_past==1)=0;


%%%Reduces insolation due to atmospheric thickness
%Nanoparticle properites
nanoparticlemass=mass;
nanoparticlecolumndensity=nanoparticlemass/(sum(sum(INTERPSTACK.gridarea*1e10)));
rayleighscatter=(550)^-4;
rayleighabsorb=(550)^-1;
rayleightextinction=rayleighscatter+rayleighabsorb;
nanoparticleopticaldepth=(nanoparticlecolumndensity*rayleightextinction)/(4*2.7*10e-9);

%Summed depths

%Reduces biomass due to asteroid impact & effects
culledmaterial(land_past==1)=1;
culledmaterial(truefires==1 & struck==0)=0.5;
culledmaterial(struck==1 & truefires==1)=0;
culledmaterial(land_past==0)=0;
for i = 1 : 40
    for j = 1 : 48
        if culledmaterial(i,j)~=0
            tocull=randi(200);
            if tocull>100
                tocull=tocull-100;
            else
                tocull=-tocull;
            end
            tocull=tocull/1000;
            culledmaterial(i,j)=culledmaterial(i,j)+tocull;
        end
    end
end
culledmaterial(culledmaterial<0)=0;
totalprevbiomass=(sum( sum( toburn .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' )));
%plantgridarea(land_past==0)=0;
totalplantgridarea=sum(sum(plantgridarea*1e10));
burntmaterial=(toburn.*(1-culledmaterial));
biomasstocarry=toburn.*culledmaterial;
charcoalmaterial=(sum( sum( biomasstocarry .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' )))*(1/3);
Crelease=(sum( sum( burntmaterial .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' ))/12);
aciniformcarbonquant=Crelease*12*0.03;
carbonvals=[Crelease*0.43 aciniformcarbonquant charcoalmaterial];
writematrix(carbonvals,"carbonoutput.csv")
culledmaterial(culledmaterial~=0)=0;
%culledmaterial(topologymap>1e3)=0.001;
culledmaterial(runoffmap>1.2e3)=0.001;
growthmap(culledmaterial~=0)=SMITE_asteroidparams(1);
%Maps wildfires
firemap=potentialfires ;
firemap(potentialfires==1 & land_past==1)=3;
firemap(potentialfires==0 & land_past==0)=0;
firemap(potentialfires==0 & land_past==1)=1;
firemap(potentialfires==1 & land_past==0)=2;
firemap(potentialfires==1 & struck==1)=4;

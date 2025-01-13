%%%%%%%%%%%%%%%%SMITE

global biopars
global culledmaterial
global INTERPSTACK
global toburn
global Crelease
global biomasstocarry
global growthmap
global SMITE_asteroidparams

GRID_AREA_km2=INTERPSTACK.gridarea;
key_time = max ( INTERPSTACK.time( (INTERPSTACK.time - SMITE_asteroidparams(1)) <= 0 ) ) ;
key_index = find( INTERPSTACK.time == key_time ) ;
land_past = INTERPSTACK.land(:,:,key_index) ; 
potentialfires = INTERPSTACK.wildfires(:,:) ;
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
        toaddy=[struckcellsy(j)-1,struckcellsy(j),struckcellsy(j)+1];
    end
    if struckcellsx(j)==29
        toaddx=[28,29,27,26];
    elseif struckcellsx(j)==28
        toaddx=[26,27,28,29];
    else
        toaddx=[struckcellsx(j)-1,struckcellsx(j),struckcellsx(j)+1];
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
culledmaterial(truefires==1 & struck==0)=0.03;
culledmaterial(struck==1 & truefires==1)=0;
culledmaterial(land_past==0)=0;
for i = 1 : 40
    for j = 1 : 48
        if culledmaterial(i,j)~=0
            culledmaterial(i,j)=culledmaterial(i,j)-(randi(200)/1e3);
        end
    end
end
culledmaterial(culledmaterial<0)=0;
growthmap(culledmaterial~=0)=SMITE_asteroidparams(1);
burntmaterial=(toburn.*(1-culledmaterial));
biomasstocarry=toburn.*culledmaterial;
Crelease=(sum( sum( burntmaterial .* ( GRID_AREA_km2 * 1e6 ), 'omitnan' ))/12);
%Maps wildfires
firemap=potentialfires ;
firemap(potentialfires==1 & land_past==1)=3;
firemap(potentialfires==0 & land_past==0)=0;
firemap(potentialfires==0 & land_past==1)=1;
firemap(potentialfires==1 & land_past==0)=2;
firemap(potentialfires==1 & struck==1)=4;

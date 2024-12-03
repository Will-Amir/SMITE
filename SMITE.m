%%%%%%%%%%%%%%%%SMITE

global asteroidtimes
global biopars
global culledmaterial
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
    asteroidtimes=valstorun(1);
    lats=valstorun(2);
    longs=valstorun(3);
    powers=valstorun(4);
    for i = 1:length(asteroidtimes)
        struck(lats(i),longs(i))=1;
    end
end
struck=circshift(struck, [0 29]);
%Generates surrounding cells for wildires
struckcells = find(struck==1) ;
struckcellsx = floor(struckcells/biopars.x_lon) ;
struckcellsy = rem(struckcells,biopars.x_lon) ;
for j = 1:length(asteroidtimes)
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
    if powers(j)>=98
        potentialfires(:,:) = 1;
    else
        potentialfires(toaddy,toaddx) = 1;
    end
end
%Checks that wildfires are on land
truefires=potentialfires ;
truefires(potentialfires==1 & land_past==1)=1 ;
truefires(potentialfires==0 & land_past==1)=0;

%Reduces biomass due to asteroid impact & effects
culledmaterial(land_past==1)=1;
culledmaterial(truefires==1 & struck==0)=0.4;
culledmaterial(struck==1 & truefires==1)=0.2;
culledmaterial(land_past==0)=0;

%Maps wildfires
firemap=potentialfires ;
firemap(potentialfires==1 & land_past==1)=3;
firemap(potentialfires==0 & land_past==0)=0;
firemap(potentialfires==0 & land_past==1)=1;
firemap(potentialfires==1 & land_past==0)=2;

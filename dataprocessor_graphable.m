global SMITE_asteroidparams
clear all
SMITE_asteroidparams(1)=-66.4;
%%% Set up arrays
biomassquant=[];
co2quant=[];
tempquant=[];
times=[];


%%%% Read data
runnums=50;
timetoextend=50000;
hostdir='C:\Users\py20waa\OneDrive - University of Leeds\Masters\Outputs\MATLAB runs\SMITE paper runs\100 year regrowth time';
for x = 1:runnums
    newdir=append("\Run",string(x));
    changedir=append(hostdir,newdir);
    cd (changedir)
    parameters=readtable("runoutput.csv");
    timeindex = find((parameters.time/1e6)==(SMITE_asteroidparams(1)))-1;

    %%Times
    timestoconsider=parameters.time(timeindex:(timeindex+timetoextend));
    newtimes=cat(2,times,timestoconsider);
    times=newtimes;

    %%Biomasses
    biomasstoconsider=parameters.biomass_tot(timeindex:(timeindex+timetoextend));
    newbiomass=cat(2,biomassquant,biomasstoconsider);
    biomassquant=newbiomass;

    %%CO2
    co2toconsider=(parameters.RCO2(timeindex:(timeindex+timetoextend))).*280;
    newco2=cat(2,co2quant,co2toconsider);
    co2quant=newco2;

    %%Temp
    temptoconsider=parameters.tempC(timeindex:(timeindex+timetoextend));
    newtemp=cat(2,tempquant,temptoconsider);
    tempquant=newtemp;
end

%%%Calculate data
%%Biomass
biomassmeans=mean(biomassquant,2);
biomassmedians=median(biomassquant,2);
biomassstdevs=std(biomassquant,0,2);

%%CO2
CO2means=mean(co2quant,2);
CO2medians=median(co2quant,2);
CO2stdevs=std(co2quant,0,2);

%%Temp
tempmeans=mean(tempquant,2);
tempmedians=median(tempquant,2);
tempstddevs=std(tempquant,0,2);

%%%Output data
%Collate data
matrixheaders=["Biomass: mean" "Median" "STDDEV" "CO2: mean" "Median" "STDDEV" "Temp: mean" "Median" "STDDEV" "Time"];
matrixdata=[biomassmeans biomassmedians biomassstdevs CO2means CO2medians CO2stdevs tempmeans tempmedians tempstddevs times(:,1)];
matrixtowrite=vertcat(matrixheaders,matrixdata);
cd (hostdir)
writematrix(matrixtowrite,"timedependentdata.csv")
cd("C:\Users\py20waa\Documents\MATLAB\SCI-FId-SMITE")
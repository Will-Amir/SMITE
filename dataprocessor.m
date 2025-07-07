global SMITE_asteroidparams


%%%% Read data
runnums=50;
hostdir='C:\Users\py20waa\OneDrive - University of Leeds\Masters\Outputs\MATLAB runs\SMITE paper runs\100 year regrowth time';
recoverytime=[];
preimpactbiomass=[];
peakbiomass=[];
co2releasequant=[];
aciniformquant=[];
charcoalquant=[];
peakco2=[];
peakco2time=[];
peaktemp=[];
peaktemptime=[];
for x = 1:runnums
    newdir=append("\Run",string(x));
    changedir=append(hostdir,newdir);
    cd (changedir)
    %%%Carbon value reader
    carbonvals=readmatrix("carbonoutput.csv");
    co2releasequant(end+1)=carbonvals(1);
    aciniformquant(end+1)=carbonvals(2);
    charcoalquant(end+1)=carbonvals(3);
    %%%Recovery times and peak value reader
    %{
    parameters=readtable("runoutput.csv");
    timeindex = find((parameters.time/1e6)==(SMITE_asteroidparams(1)))-1;
    maximums = [max(parameters.temperature) max(parameters.RCO2)];
    maximumindex = [find(parameters.temperature==maximums(1)) find(parameters.RCO2==maximums(2))];
    percentilescale = [parameters.temperature(maximumindex(1)) parameters.RCO2(maximumindex(2))*280 parameters.biomass_tot(timeindex)];
    
    maxCO2timeindex=parameters.time(maximumindex(2));
    
    maxtemptimeindex=parameters.time(maximumindex(1));
    
    biomasspercentilescale=(parameters.biomass_tot(:)./percentilescale(3))*100;
    possiblebiomasses=find(biomasspercentilescale>=100);
    recoverytimeindiv=min(possiblebiomasses(possiblebiomasses>timeindex));
    peakbiomassindiv=max(possiblebiomasses(possiblebiomasses>timeindex));


    recoverytime(end+1)=recoverytimeindiv;
    peakco2(end+1)=maximums(2)*280;
    peakco2time(end+1)=maximumindex(2);
    peaktemp(end+1)=maximums(1);
    peaktemptime(end+1)=maximumindex(1);
    preimpactbiomass(end+1)=parameters.biomass_tot(timeindex);
    peakbiomass(end+1)=peakbiomassindiv;
end
}

%%%Process data
%Recovery time
recoverytimemean=mean(recoverytime);
recoverytimedev=std(recoverytime);
recoverytimerange=range(recoverytime);
recoverytimemedian=median(recoverytime);
%Preimpact biomass
preimpactbiomassmean=mean(preimpactbiomass);
preimpactbiomassdev=std(preimpactbiomass);
preimpactbiomassrange=range(preimpactbiomass);
preimpactbiomassmedian=median(preimpactbiomass);
%Peak biomass
peakbiomassmean=mean(peakbiomass);
peakbiomassdev=std(peakbiomass);
peakbiomassrange=range(peakbiomass);
peakbiomassmedian=median(peakbiomass);
%CO2 release
co2releasequantmean=mean(co2releasequant);
co2releasequantdev=std(co2releasequant);
co2releasequantrange=range(co2releasequant);
co2releasequantmedian=median(co2releasequant);
%Aciniform carbon
aciniformquantmean=mean(aciniformquant);
aciniformquantdev=std(aciniformquant);
aciniformquantrange=range(aciniformquant);
aciniformquantmedian=median(aciniformquant);
%Charcoal carbon
charcoalquantmean=mean(charcoalquant);
charcoalquantdev=std(charcoalquant);
charcoalquantrange=range(charcoalquant);
charcoalquantmedian=median(charcoalquant);
%Peak CO2
peakco2mean=mean(peakco2);
peakco2dev=std(peakco2);
peakco2range=range(peakco2);
peakco2median=median(peakco2);
%Peak CO2 time
peakco2timemean=mean(peakco2time);
peakco2timedev=std(peakco2time);
peakco2timerange=range(peakco2time);
peakco2timemedian=median(peakco2time);
%Peak temp
peaktempmean=mean(peaktemp);
peaktempdev=std(peaktemp);
peaktemprange=range(peaktemp);
peaktempmedian=median(peaktemp);
%Peak temp time
peaktemptimemean=mean(peaktemptime);
peaktemptimedev=std(peaktemptime);
peaktemptimerange=range(peaktemptime);
peaktemptimemedian=median(peaktemptime);


%%%Output data
%Collate data
matrixtowrite=["Recovery time" NaN NaN NaN ; recoverytimemean recoverytimedev recoverytimerange recoverytimemedian ; "Preimpact biomass" NaN NaN NaN ; preimpactbiomassmean preimpactbiomassdev preimpactbiomassrange preimpactbiomassmedian ; "Peak biomass" NaN NaN NaN ; peakbiomassmean peakbiomassdev peakbiomassrange peakbiomassmedian ; "CO2 release" NaN NaN NaN ; co2releasequantmean co2releasequantdev co2releasequantrange co2releasequantmedian ; "Aciniform carbon" NaN NaN NaN ; aciniformquantmean aciniformquantdev aciniformquantrange aciniformquantmedian ; "Charcoal carbon" NaN NaN NaN ; charcoalquantmean charcoalquantdev charcoalquantrange charcoalquantmedian ; "Peak CO2 level" NaN NaN NaN ; peakco2mean peakco2dev peakco2range peakco2median ; "Peak CO2 time" NaN NaN NaN ; peakco2timemean peakco2timedev peakco2timerange peakco2timemedian ; "Peak temperature" NaN NaN NaN ; peaktempmean peaktempdev peaktemprange peaktempmedian ; "Peak temperature time" NaN NaN NaN ; peaktemptimemean peaktemptimedev peaktemptimerange peaktemptimemedian];
cd (hostdir)
writematrix(matrixtowrite,"processeddata.csv")
cd("C:\Users\py20waa\Documents\MATLAB\SCI-FId-SMITE")
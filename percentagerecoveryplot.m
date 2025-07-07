%%%% Percentage recovery plots
global state
global SMITE_asteroidparams
parameters=state;%readtable("runoutput.csv");
timeindex = find(parameters.time==(SMITE_asteroidparams(1)*1e6))-1;
maximums = [max(parameters.temperature) max(parameters.RCO2)];
maximumindex = [find(parameters.temperature==max(parameters.temperature)) find(parameters.RCO2==max(parameters.RCO2))];
percentilescale = [parameters.temperature(maximumindex(1)) parameters.RCO2(maximumindex(2))*280 parameters.biomass_tot(timeindex)];

CO2percentile=percentilescale(2)-parameters.RCO2(timeindex)*280;
CO2percentilescale=((percentilescale(2)-parameters.RCO2(:)*280)./CO2percentile)*100;
maxCO2timeindex=parameters.time(maximumindex(2));

temppercentile=percentilescale(1)-parameters.temperature(timeindex);
temppercentilescale=((percentilescale(1)-parameters.temperature(:))./temppercentile)*100;
maxtemptimeindex=parameters.time(maximumindex(1));

biomasspercentilescale=(parameters.biomass_tot(:)./percentilescale(3))*100;
figure()
subplot(1,3,1)
hold on
xlim([maxCO2timeindex -65.4e6])
xlabel('Time (Ma)')
ylim([0 100])
ylabel('CO2 recovery percentile (%)')
title("CO2 recovery")
yline(10,"r--")
yline(20,"b:")
yline(50,"m-.")
plot(parameters.time,CO2percentilescale)

subplot(1,3,2)
hold on
xlim([maxtemptimeindex -65.4e6])
xlabel('Time (Ma)')
ylim([0 100])
ylabel('Temperature recovery percentile (%)')
title("Temperature recovery")
yline(10,"r--")
yline(20,"b:")
yline(50,"m-.")
plot(parameters.time,temppercentilescale)

subplot(1,3,3)
hold on
xlim([SMITE_asteroidparams(1)*1e6 -65.4e6])
xlabel('Time (Ma)')
ylim([0 100])
ylabel('Biomass recovery percentile (%)')
title("Biomass recovery")
yline(10,"r--")
yline(20,"b:")
yline(50,"m-.")
plot(parameters.time,biomasspercentilescale)

savefig("recoveryplots")
clear all
load('data/updated_historical_isotopes.mat')
cd("C:\Users\py20waa\Documents\MATLAB\SCI-FId-SMITE\DataToWorkon\")



parameters100=readtable("timedependentdata_100yr.csv");
parameters1000=readtable("timedependentdata_1000yr.csv");
parameters10000=readtable("timedependentdata_10000yr.csv");
parameters100000=readtable("timedependentdata_100000yr.csv");

x1 = parameters100.Time';% + 66.4e6;
y1 = parameters100.Temp_Mean';
stddev1 = parameters100.STDDEV_2';

uppery1 = (y1+2*stddev1);
lowery1 = (y1-2*stddev1);

xdevs1 = [x1, fliplr(x1)];
ydevs1 = [uppery1, fliplr(lowery1)];

x2 = parameters1000.Time';% + 66.4e6;
y2 = parameters1000.Temp_Mean';
stddev2 = parameters1000.STDDEV_2';

uppery2 = (y2+2*stddev2);
lowery2 = (y2-2*stddev2);

xdevs2 = [x2, fliplr(x2)];
ydevs2 = [uppery2, fliplr(lowery2)];

x3 = parameters10000.Time';
y3 = parameters10000.Temp_Mean';
stddev3 = parameters10000.STDDEV_2';

uppery3 = (y3+2*stddev3);
lowery3 = (y3-2*stddev3);

xdevs3 = [x3, fliplr(x3)];
ydevs3 = [uppery3, fliplr(lowery3)];

x4 = parameters100000.Time';
y4 = parameters100000.Temp_Mean';
stddev4 = parameters100000.STDDEV_2';

uppery4 = (y4+2*stddev4);
lowery4 = (y4-2*stddev4);

xdevs4 = [x4, fliplr(x4)];
ydevs4 = [uppery4, fliplr(lowery4)];

timeconvertedoxygen=(d18o_x_highfid*1e6);
firsttimepoint=find(timeconvertedoxygen<-66.4e6,1);
temptonormalise=avgsurftemps_highfid(firsttimepoint);
toremove=temptonormalise-y1(1);
normalisedtemps=avgsurftemps_highfid-toremove+1;



figure('color','white')
hold on 
xlim([-66.4e6 -66.35e6]);

plot(x1,y1,'color',[(220/256) (5/256) (12/256)],'LineWidth',2)
plot(x2,y2,'color',[(25/256) (101/256) (176/256)],'LineWidth',2)
%plot(x3,y3,'color',[(136/256) (46/256) (114/256)],'LineWidth',2)
%plot(x4,y4,'color',[(246/256) (193/256) (65/256)],'LineWidth',2)
plot((d18o_x_highfid*1e6),normalisedtemps,'.','color',[(78/256) (178/256) (101/256)],'MarkerSize',10)
p = fill(xdevs1, ydevs1, [(220/256) (5/256) (12/256)], 'FaceAlpha',0.15,'EdgeColor','none');
p = fill(xdevs2, ydevs2, [(25/256) (101/256) (176/256)], 'FaceAlpha',0.15,'EdgeColor','none');
%p = fill(xdevs3, ydevs3, [(136/256) (46/256) (114/256)], 'FaceAlpha',0.15,'EdgeColor','none');
%p = fill(xdevs4, ydevs4, [(246/256) (193/256) (65/256)], 'FaceAlpha',0.15,'EdgeColor','none');
h = gca;
h.XTickMode = 'auto';
h.XTickLabel = string((h.XTick+66.4e6));
h.XLabel.String = 'Time since K-Pg impact event (Yr)';
h.YTickMode = 'auto';
h.YLabel.String = 'Atmospheric CO_2 (ppm)';
lgd=legend('100 year regrowth time ','1000 year regrowth time', 'Normalised Historical data', '2σ value', '2σ value');
h.FontSize=20;
lgd.FontSize=20;

hold off
cd ("C:\Users\py20waa\Documents\MATLAB\SCI-FId-SMITE")


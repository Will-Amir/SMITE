filetoread = "aba6853_tables_s8_s34.xlsx";
excelfile = readtable(filetoread,"Sheet","Table S33");
d13c_x_highfid = -(transpose(excelfile.age_tuned));
d13c_y_highfid = transpose(excelfile.benthicD13CVPDB);
d18o_x_highfid = -(transpose(excelfile.age_tuned));
d18o_y_highfid = transpose(excelfile.benthicD18OVPDB);
save("updated_historical_isotopes.mat","d13c_x_highfid","d13c_y_highfid","d18o_x_highfid","d18o_y_highfid");

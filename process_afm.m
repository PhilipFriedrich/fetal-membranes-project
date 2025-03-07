clc, clear, close all

% define base folder
base_folder = '';

% load the compiled data
load([base_folder 'afm_data.mat']);

% get indices for different sample types
amni = fm_data.('layer (1-a/0-c)');
chori = ~amni;
normal = fm_data.('delivery (1-n/0-ps)');
sectio = ~normal;

%%

% define plot style
paint = [212,175,55]./255;
alpha = 0.8;
markeralpha = 0.6;
titlename = 'AFM';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';


for j = 1:10
    switch j
        case 1
            % plot stiffness
            data1 = fm_data.YM_log(amni & normal);
            data2 = fm_data.YM_log(amni & sectio);
            labelname = 'Young''s modulus (Pa)';
            name = 'YM_amnion';
        case 2
            % plot storage modulus
            data1 = fm_data.G1_log(amni & normal) ./ 10^3;
            data2 = fm_data.G1_log(amni & sectio) ./ 10^3;
            labelname = 'Shear storage modulus G'' (kPa)';
            name = 'G1_amnion';
        case 3
            % plot loss modulus
            data1 = fm_data.G2_log(amni & normal) ./ 10^3;
            data2 = fm_data.G2_log(amni & sectio) ./ 10^3;
            labelname = 'Shear loss modulus G'''' (kPa)';
            name = 'G2_amnion';
        case 4
            % plot resistance
            data1 = fm_data.Gs(amni & normal) ./ 10^3;
            data2 = fm_data.Gs(amni & sectio) ./ 10^3;
            labelname = 'Resistance |G*| (kPa)';
            name = 'Gs_amnion';
        case 5
            % plot fluidity
            data1 = fm_data.phi(amni & normal);
            data2 = fm_data.phi(amni & sectio);
            labelname = 'Loss factor tan(\delta)';
            name = 'lossfactor_amnion';
        case 6
            % plot stiffness
            data1 = fm_data.YM_log(chori & normal);
            data2 = fm_data.YM_log(chori & sectio);
            labelname = 'Young''s modulus (Pa)';
            name = 'YM_chorion';
        case 7
            % plot storage modulus
            data1 = fm_data.G1_log(chori & normal) ./ 10^3;
            data2 = fm_data.G1_log(chori & sectio) ./ 10^3;
            labelname = 'Shear storage modulus G'' (kPa)';
            name = 'G1_chorion';
        case 8
            % plot loss modulus
            data1 = fm_data.G2_log(chori & normal) ./ 10^3;
            data2 = fm_data.G2_log(chori & sectio) ./ 10^3;
            labelname = 'Shear loss modulus G'''' (kPa)';
            name = 'G2_chorion';
        case 9
            % plot resistance
            data1 = fm_data.Gs(chori & normal) ./ 10^3;
            data2 = fm_data.Gs(chori & sectio) ./ 10^3;
            labelname = 'Resistance |G*| (kPa)';
            name = 'Gs_chorion';
        case 10
            % plot fluidity
            data1 = fm_data.phi(chori & normal);
            data2 = fm_data.phi(chori & sectio);
            labelname = 'Loss factor tan(\delta)';
            name = 'lossfactor_chorion';
    end

    disp(['\%%%%% testing ' name ' for normal distribution %%%%%%'])
    normalflag = [false, false];
    [SW, ~] = swft(data1);
    if SW{2,7} > 0.05
        disp(['data 1 normally distributed'])
        normalflag(1) = true;
    else
        disp(['data 1 not normally distributed'])
    end  
    [SW, ~] = swft(data2);
    if SW{2,7} > 0.05
        disp(['data 2 normally distributed'])
        normalflag(2) = true;
    else
        disp(['data 2 not normally distributed'])
    end

    [fig, p, sigtest, median_values, n] = make_boxplot_fetal_membranes(data1, data2, titlename, labelname, paint, alpha, markeralpha, colorname1, colorname2, groupname1, groupname2);
    
    exportgraphics(gcf,'', 'Resolution', 600)
    
    StatTable = create_statistics_table_fetal_membranes(data1, data2, n, p, median_values, sigtest);
    writetable(StatTable,'','Delimiter','\t') 
    
    close
end

%%
% now plot amnion adn chorion in one diagram

% define plot style
paint = [212,175,55]./255;
alpha = 0.8;
markeralpha = 0.6;
titlename = 'AFM';
groupname1 = 'Amnion';
groupname2 = 'Chorion';
colorname1 = 'SD';
colorname2 = 'CS';


for j = 1:5
    switch j
        case 1
            % plot stiffness
            data1 = fm_data.YM_log(amni & normal);
            data2 = fm_data.YM_log(amni & sectio);
            data3 = fm_data.YM_log(chori & normal);
            data4 = fm_data.YM_log(chori & sectio);
            labelname = 'Young''s modulus (Pa)';
            name = 'YM';
        case 2
            % plot storage modulus
            data1 = fm_data.G1_log(amni & normal) ./ 10^3;
            data2 = fm_data.G1_log(amni & sectio) ./ 10^3;
            data3 = fm_data.G1_log(chori & normal) ./ 10^3;
            data4 = fm_data.G1_log(chori & sectio) ./ 10^3;
            labelname = 'Shear storage modulus G'' (kPa)';
            name = 'G1';
        case 3
            % plot loss modulus
            data1 = fm_data.G2_log(amni & normal) ./ 10^3;
            data2 = fm_data.G2_log(amni & sectio) ./ 10^3;
            data3 = fm_data.G2_log(chori & normal) ./ 10^3;
            data4 = fm_data.G2_log(chori & sectio) ./ 10^3;
            labelname = 'Shear loss modulus G'''' (kPa)';
            name = 'G2';
        case 4
            % plot resistance
            data1 = fm_data.Gs(amni & normal) ./ 10^3;
            data2 = fm_data.Gs(amni & sectio) ./ 10^3;
            data3 = fm_data.Gs(chori & normal) ./ 10^3;
            data4 = fm_data.Gs(chori & sectio) ./ 10^3;
            labelname = 'Resistance |G*| (kPa)';
            name = 'Gs';
        case 5
            % plot fluidity
            data1 = fm_data.phi(amni & normal);
            data2 = fm_data.phi(amni & sectio);
            data3 = fm_data.phi(chori & normal);
            data4 = fm_data.phi(chori & sectio);
            labelname = 'Loss factor tan(\delta)';
            name = 'lossfactor';
    end

[fig, p, sigtest, median_values, n, nature_data] = make_boxplot_fetal_membranes_v2(data1, data2, data3, data4, titlename, labelname, paint, alpha, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)

StatTable = create_statistics_table_fetal_membranes_v2(data1, data2, data3, data4, n, p, median_values, sigtest);
writetable(StatTable,'','Delimiter','\t') 

save('', 'nature_data')

close
end

disp('done.')


% load patient information
patienttable = readtable('fetal_membranes_patient_info.xlsx');

idx = 0;
for i = 1:length(fm_data.number)
    try
        idx(i,1) = find(patienttable.ProbenNr == fm_data.number(i));
    catch
        disp(['>>> sample ' num2str(fm_data.number(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% append fm_data with patient data
fm_data.SST = patienttable.SST(idx);
fm_data.KU = patienttable.KU_cm_(idx);
fm_data.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
fm_data.Sex = patienttable.Geschlecht_1_m_(idx);
fm_data.Gw_Kind = patienttable.GewichtKind_g_(idx);
fm_data.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
fm_data.Gr_Mutter = patienttable.Gr__e_cm_(idx);
fm_data.BMI = patienttable.BMIVorSS(idx);
fm_data.Nikotin = patienttable.Nikotin(idx);


%%

% add ssw plot 
% define plot style
paint = [212,175,55]./255;
alpha = 0.8;
markeralpha = 0.8;
titlename = 'AFM';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';

% SSW
for j = 1:2
    switch j
        case 1
            % plot resistance
            data1 = fm_data.Gs(amni & normal) ./ 10^3;
            data2 = fm_data.Gs(amni & sectio) ./ 10^3;
            time1 = fm_data.SST(amni & normal);
            time2 = fm_data.SST(amni & sectio);
            labelname = 'Resistance |G*| (kPa)';
            name = 'ssw_amnion';
        case 2
            % plot resistance
            data1 = fm_data.Gs(chori & normal) ./ 10^3;
            data2 = fm_data.Gs(chori & sectio) ./ 10^3;
            time1 = fm_data.SST(chori & normal);
            time2 = fm_data.SST(chori & sectio);
            labelname = 'Resistance |G*| (kPa)';
            name = 'ssw_chorion';
    end

 [fig, n] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)
close
end
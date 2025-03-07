clc, clear, close all

% define base folder
base_folder = '';

% load the compiled data
load([base_folder 'rheo_data.mat']);

rheo_data = dataRheoCorrected;

% get idneices for frequency of interest
freq_idx = (rheo_data.FS_Frequency{1} == 1);

%%

% get data of interest 
for i = 1:height(rheo_data)
    G1(i,1) = rheo_data.FS_StorageModulus{i}(freq_idx);
    G2(i,1) = rheo_data.FS_LossModulus{i}(freq_idx);
    Gs(i,1) = rheo_data.FS_ComplexShearModulus{i}(freq_idx);
    phi(i,1) = rheo_data.FS_LossFactor{i}(freq_idx);
    ssw(i,1) = str2double(rheo_data.ssw{i}) * 7 + str2double(rheo_data.sst{i});
    number(i,1) = str2double(rheo_data.probeNr{i});
end

% remove G1 too high (53, 41, 47, 48, 42) or too low (29) as they need
% to be excluded from analysis
del_idx = (number == [53, 41, 47, 48, 42, 29]);
del_idx = ~sum(del_idx,2);

G1 = G1(del_idx);
G2 = G2(del_idx);
Gs = Gs(del_idx);
phi = phi(del_idx);
ssw = ssw(del_idx);
number = number(del_idx);

% get type of sample
normal = (rheo_data.fall(del_idx) == 'spontan');
sectio = (rheo_data.fall(del_idx) == 'sectio');

% define plot style
paint = [39,107,97]./255;
alpha = 0.8;
markeralpha = 0.6;
titlename = 'Rheometer';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';


for j = 1:4
    switch j
        case 1
            % plot storage modulus
            data1 = G1(normal) ./ 10^3;
            data2 = G1(sectio) ./ 10^3;
            labelname = 'Shear storage modulus G'' (kPa)';
            name = 'G1';
        case 2
            % plot loss modulus
            data1 = G2(normal) ./ 10^3;
            data2 = G2(sectio) ./ 10^3;
            labelname = 'Shear loss modulus G'''' (kPa)';
            name = 'G2';
        case 3
            % plot resistance
            data1 = Gs(normal) ./ 10^3;
            data2 = Gs(sectio) ./ 10^3;
            labelname = 'Resistance |G*| (kPa)';
            name = 'Gs';
        case 4
            % plot fluidity
            data1 = phi(normal);
            data2 = phi(sectio);
            labelname = 'Loss factor tan(\delta)';
            name = 'lossfactor';
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


[fig, p, sigtest, median_values, n, nature_data] = make_boxplot_fetal_membranes(data1, data2, titlename, labelname, paint, alpha, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)

StatTable = create_statistics_table_fetal_membranes(data1, data2, n, p, median_values, sigtest);
writetable(StatTable,'','Delimiter','\t') 

save('', 'nature_data')

close
end


%%

% load patient information
patienttable = readtable('fetal_membranes_patient_info.xlsx');

rheo_data.number = str2double(rheo_data.probeNr);
idx = 0;
for i = 1:length(rheo_data.number)
    try
        idx(i,1) = find(patienttable.ProbenNr == rheo_data.number(i));
    catch
        disp(['>>> sample ' num2str(rheo_data.number(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% append rheo_data with patient data
rheo_data.SST = patienttable.SST(idx);
rheo_data.KU2 = patienttable.KU_cm_(idx);
rheo_data.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
rheo_data.Sex = patienttable.Geschlecht_1_m_(idx);
rheo_data.Gw_Kind = patienttable.GewichtKind_g_(idx);
rheo_data.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
rheo_data.Gr_Mutter = patienttable.Gr__e_cm_(idx);
rheo_data.BMI = patienttable.BMIVorSS(idx);
rheo_data.Nikotin = patienttable.Nikotin(idx);

%%

% define plot style
paint = [39,107,97]./255;
alpha = 0.8;
markeralpha = 0.8;
titlename = 'Rheometer';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';

% plot resistance
data1 = Gs(normal) ./ 10^3;
data2 = Gs(sectio) ./ 10^3;
time1 = rheo_data.SST(normal);
time2 = rheo_data.SST(sectio);
labelname = 'Resistance |G*| (kPa)';
name = 'ssw';

[fig, n, nature_data] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)

save('', 'nature_data')


close



clc, clear, close all

% define base folder

base_folder = '';

load([base_folder 'mre_data.mat']);

mre_data = GroupList;


%%

% define plot style
paint = [204,85,0]./255;
alpha = 0.8;
markeralpha = 0.6;
titlename = 'MRE';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';


for j = 1:5
    switch j
        case 1
            % plot storage modulus
            data1 = mre_data(1).Gp1000 ./ 10^3;
            data2 = mre_data(2).Gp1000 ./ 10^3;
            labelname = 'Shear storage modulus G'' (kPa)';
            name = 'G1';
        case 2
            % plot loss modulus
            data1 = mre_data(1).Gdp1000 ./ 10^3;
            data2 = mre_data(2).Gdp1000 ./ 10^3;
            labelname = 'Shear loss modulus G'''' (kPa)';
            name = 'G2';
        case 3
            % plot resistance
            data1 = mre_data(1).Gstar1000 ./ 10^3;
            data2 = mre_data(2).Gstar1000 ./ 10^3;
            labelname = 'Resistance |G*| (kPa)';
            name = 'Gs';
        case 4
            % plot fluidity
            data1 = mre_data(1).phi;
            data2 = mre_data(2).phi;
            labelname = 'Loss factor tan(\delta)';
            name = 'losfactor';
        case 5
            % plot ADC
            data1 = mre_data(1).ADC .* 10^6; % convert to mm² instead of m²
            data2 = mre_data(2).ADC .* 10^6;
            labelname = 'Water diffusion (mm²/s)';
            name = 'ADC';
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
    
    
    % change y axis interval for loss factor only
    if j == 4
        set(gca, 'ytick', 0:0.2:3);
    end
    
    exportgraphics(gcf,'', 'Resolution', 600)
    
    StatTable = create_statistics_table_fetal_membranes(data1, data2, n, p, median_values, sigtest);
    writetable(StatTable,'','Delimiter','\t') 
    
    save('', 'nature_data')

    close
end


%%


% load patient information
patienttable = readtable('fetal_membranes_patient_info.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normal mre
idx = 0;
for i = 1:length(mre_data(1).MREnumber)
    try
        idx(i,1) = find(patienttable.ProbenNr == mre_data(1).MREnumber(i));
    catch
        disp(['>>> sample ' num2str(mre_data(1).MREnumber(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% get patient data
temp = [];
temp = array2table(patienttable.SST(idx),'VariableNames',{'SST'});
temp.KU = patienttable.KU_cm_(idx);
temp.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
temp.Sex = patienttable.Geschlecht_1_m_(idx);
temp.Gw_Kind = patienttable.GewichtKind_g_(idx);
temp.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
temp.Gr_Mutter = patienttable.Gr__e_cm_(idx);
temp.BMI = patienttable.BMIVorSS(idx);
temp.Nikotin = patienttable.Nikotin(idx);

mre_normal_table = [array2table([mre_data(1).Gp1000, mre_data(1).Gdp1000, mre_data(1).Gstar1000, mre_data(1).phi],'VariableNames',{'G1','G2','Gs','phi'}), temp];
mre_normal_sst = temp.SST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sectio mre
idx = 0;
for i = 1:length(mre_data(2).MREnumber)
    try
        idx(i,1) = find(patienttable.ProbenNr == mre_data(2).MREnumber(i));
    catch
        disp(['>>> sample ' num2str(mre_data(2).MREnumber(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% get patient data
temp = [];
temp = array2table(patienttable.SST(idx),'VariableNames',{'SST'});
temp.KU = patienttable.KU_cm_(idx);
temp.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
temp.Sex = patienttable.Geschlecht_1_m_(idx);
temp.Gw_Kind = patienttable.GewichtKind_g_(idx);
temp.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
temp.Gr_Mutter = patienttable.Gr__e_cm_(idx);
temp.BMI = patienttable.BMIVorSS(idx);
temp.Nikotin = patienttable.Nikotin(idx);

mre_sectio_table = [array2table([mre_data(2).Gp1000, mre_data(2).Gdp1000, mre_data(2).Gstar1000, mre_data(2).phi],'VariableNames',{'G1','G2','Gs','phi'}), temp];
mre_sectio_sst = temp.SST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normal adc
idx = 0;
for i = 1:length(mre_data(1).ADCnumber)
    try
        idx(i,1) = find(patienttable.ProbenNr == mre_data(1).ADCnumber(i));
    catch
        disp(['>>> sample ' num2str(mre_data(1).ADCnumber(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% get patient data
temp = [];
temp = array2table(patienttable.SST(idx),'VariableNames',{'SST'});
temp.KU = patienttable.KU_cm_(idx);
temp.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
temp.Sex = patienttable.Geschlecht_1_m_(idx);
temp.Gw_Kind = patienttable.GewichtKind_g_(idx);
temp.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
temp.Gr_Mutter = patienttable.Gr__e_cm_(idx);
temp.BMI = patienttable.BMIVorSS(idx);
temp.Nikotin = patienttable.Nikotin(idx);

adc_normal_table = [array2table([mre_data(1).ADC],'VariableNames',{'ADC'}), temp];
adc_normal_sst = temp.SST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sectio adc
idx = 0;
for i = 1:length(mre_data(2).ADCnumber)
    try
        idx(i,1) = find(patienttable.ProbenNr == mre_data(2).ADCnumber(i));
    catch
        disp(['>>> sample ' num2str(mre_data(2).ADCnumber(i)) ' not found']);
        idx(i,1) = nan;
    end
end

% get patient data
temp = [];
temp = array2table(patienttable.SST(idx),'VariableNames',{'SST'});
temp.KU = patienttable.KU_cm_(idx);
temp.Gr_Kind = patienttable.Gr__eKind_cm_(idx);
temp.Sex = patienttable.Geschlecht_1_m_(idx);
temp.Gw_Kind = patienttable.GewichtKind_g_(idx);
temp.Gw_Mu = patienttable.GewichtVorSSInKg(idx);
temp.Gr_Mutter = patienttable.Gr__e_cm_(idx);
temp.BMI = patienttable.BMIVorSS(idx);
temp.Nikotin = patienttable.Nikotin(idx);

adc_sectio_table = [array2table([mre_data(2).ADC],'VariableNames',{'ADC'}), temp];
adc_sectio_sst = temp.SST;

%%

% define plot style
paint = [204,85,0]./255;
alpha = 0.8;
markeralpha = 0.8;
titlename = 'MRE';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';

% plot resistance
data1 = mre_data(1).Gstar1000 ./ 10^3;
data2 = mre_data(2).Gstar1000 ./ 10^3;
time1 = mre_normal_sst;
time2 = mre_sectio_sst;
data1 = data1(~isnan(time1));
data2 = data2(~isnan(time2));
time1 = time1(~isnan(time1));
time2 = time2(~isnan(time2));
labelname = 'Resistance |G*| (kPa)';
name = 'ssw_Gs';
[fig, n, nature_data] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)

save('', 'nature_data')


close


% plot G1
data1 = mre_data(1).Gp1000 ./ 10^3;
data2 = mre_data(2).Gp1000 ./ 10^3;
time1 = mre_normal_sst;
time2 = mre_sectio_sst;
data1 = data1(~isnan(time1));
data2 = data2(~isnan(time2));
time1 = time1(~isnan(time1));
time2 = time2(~isnan(time2));
labelname = 'Shear storage modulus G'' (kPa)';
name = 'ssw_G1';
[fig, n] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)
close


% plot G2
data1 = mre_data(1).Gdp1000 ./ 10^3;
data2 = mre_data(2).Gdp1000 ./ 10^3;
time1 = mre_normal_sst;
time2 = mre_sectio_sst;
data1 = data1(~isnan(time1));
data2 = data2(~isnan(time2));
time1 = time1(~isnan(time1));
time2 = time2(~isnan(time2));
labelname = 'Shear loss modulus G'''' (kPa)';
name = 'ssw_G2';
[fig, n] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)
close

% plot phi
data1 = mre_data(1).phi;
data2 = mre_data(2).phi;
time1 = mre_normal_sst;
time2 = mre_sectio_sst;
data1 = data1(~isnan(time1));
data2 = data2(~isnan(time2));
time1 = time1(~isnan(time1));
time2 = time2(~isnan(time2));
labelname = 'Loss factor tan(\delta)';
name = 'ssw_lossfactor';
[fig, n] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)
close


% plot adc
data1 = mre_data(1).ADC .* 10^6;
data2 = mre_data(2).ADC .* 10^6;
time1 = adc_normal_sst;
time2 = adc_sectio_sst;
data1 = data1(~isnan(time1));
data2 = data2(~isnan(time2));
time1 = time1(~isnan(time1));
time2 = time2(~isnan(time2));
labelname = 'Water diffusion (mm²/s)';
name = 'ssw_ADC';
[fig, n] = make_scatter_fetal_membranes(data1, data2, time1, time2, titlename, labelname, paint, markeralpha, colorname1, colorname2, groupname1, groupname2);

exportgraphics(gcf,'', 'Resolution', 600)
close



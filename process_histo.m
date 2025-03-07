clc, clear, close all

% define base folder
base_folder = '';

load([base_folder, 'histo_data_1.mat']);

% clean_variables from bad processed immages
out_idx = [3, 7, 8, 9, 17, 28, 33, 34, 37, 38, 39, 40, 41, 44, 21, 52, 57, 4];
clear_idx = [];
for i = 1:length(saveName)
    if ismember(i, out_idx)
        clear_idx(i) = true;
    else
        clear_idx(i) = false;
    end
end
saveName = saveName(~clear_idx);
saveTable = saveTable(~clear_idx,:);

% arabisch ist normal = spontan
normal.observableNames = observableNames;
normal.saveName = saveName';
normal.saveTable = saveTable;

clear observableNames saveName saveTable out_idx clear_idx

load([base_folder, 'histo_data_2.mat']);

% clean_variables from bad processed immages
out_idx = [5, 6, 9, 11, 32, 33, 55, 9, 7, 10, 13, 14, 15];
clear_idx = [];
for i = 1:length(saveName)
    if ismember(i, out_idx)
        clear_idx(i) = true;
    else
        clear_idx(i) = false;
    end
end
saveName = saveName(~clear_idx);
saveTable = saveTable(~clear_idx,:);

% römisch ist sectio
sectio.observableNames = observableNames;
sectio.saveName = saveName';
sectio.saveTable = saveTable;
clear observableNames saveName saveTable out_idx clear_idx

%%

% define plot style
paint = [70, 130, 180]./255;
alpha = 0.8;
markeralpha = 0.6;
titlename = '';
groupname1 = '';
groupname2 = '';
colorname1 = 'SD';
colorname2 = 'CS';

% plot typical boxplots
for j = 1:3
    switch j
        case 1
            % plot fiber alignment (=Mean nematic order)
            data1 = cell2mat(normal.saveTable(:,1));
            data2 = cell2mat(sectio.saveTable(:,1));
            labelname = 'Global fiber alignment';
            name = 'Histo_order';
        case 2
            % plot thickness (=Length scale of ICT)
            data1 = cell2mat(normal.saveTable(:,3)) .* 2; % short axis starts at middle of "circle", so use factor 2 to get thickness/diameter
            data2 = cell2mat(sectio.saveTable(:,3)) .* 2; 
            labelname = 'ICT thickness (µm)';
            name = 'Histo_thickness'; 
        case 3
            % plot layer density (=ECM area fraction)
            data1 = cell2mat(normal.saveTable(:,4));
            data2 = cell2mat(sectio.saveTable(:,4));
            labelname = 'Global fiber area fraction';
            name = 'Histo_density'; 
    end

    [fig, p, sigtest, median_values, n, nature_data] = make_boxplot_fetal_membranes(data1, data2, titlename, labelname, paint, alpha, markeralpha, colorname1, colorname2, groupname1, groupname2);
    
    if j == 3 || j == 1
        % an area fraction above 1 does not exist
        axeshandles = findobj(fig, 'Type', 'axes');
        axeshandles.YTick = axeshandles.YTick(axeshandles.YTick <= 1.05);
    end

    exportgraphics(gcf,'', 'Resolution', 600)
    
    StatTable = create_statistics_table_fetal_membranes(data1, data2, n, p, median_values, sigtest);
    writetable(StatTable,'','Delimiter','\t') 
    
    save('', 'nature_data')

    close
end

%%

% extract distances and define bin edges
distances = cell2mat(normal.saveTable(1,6));
bin_edges = (0.02:0.02:1);

% normalize distances for order and density and adjust to bins
for i = 1:4
    switch i
        case 1
            % order
            datatable = cell2mat(normal.saveTable(:,7));
            ci_data = cell2mat(normal.saveTable(:,13));
        case 2
            % area fraction
            datatable = cell2mat(normal.saveTable(:,9));
            ci_data = cell2mat(normal.saveTable(:,12));
        case 3
            datatable = cell2mat(sectio.saveTable(:,7));
            ci_data = cell2mat(sectio.saveTable(:,13));
        case 4
            datatable = cell2mat(sectio.saveTable(:,9));
            ci_data = cell2mat(sectio.saveTable(:,12));
    end

    % get amount of individual image sections
    ImageNumber = height(datatable);

    % some of the denstiy values are 0 and can be set to NaN;
    datatable(datatable == 0) = NaN;
    ci_data(ci_data == 0) = NaN;

    % calculate normlaized distances
    dist = (datatable ./ datatable) .* distances;
    norm_dist = dist ./ max(dist,[],2);

    % bin the data
    binned = discretize(norm_dist, bin_edges);

    % create empty array
    nan_data = nan(size(binned));
    nan_ci = nan_data;

    % sort datatable into new array acording to binning position "binned"
    for k = 1:height(datatable)
        for kk = 1:length(datatable(k,:))
            try
                nan_data(k,binned(k,kk)) = datatable(k,kk);
                nan_ci(k,binned(k,kk)) = ci_data(k,kk);
            catch
            end
        end
    end

    % calculate the means
    data_mean = nanmean(nan_data, 1);
    ci_mean = nanmean(nan_ci, 1);

    % save data
    switch i
        case 1
            normal.mean_order_distance = data_mean;
            normal.mean_order_CI = ci_mean;
            normal.num = ImageNumber;
        case 2
            normal.mean_density_distance = data_mean;
            normal.mean_density_CI = ci_mean;
        case 3
            sectio.mean_order_distance = data_mean;
            sectio.mean_order_CI = ci_mean;
            sectio.num = ImageNumber;
        case 4
            sectio.mean_density_distance = data_mean;
            sectio.mean_density_CI = ci_mean;
    end
end


% plot distance dependencies
for j = 1:2
    switch j
        case 1
            % plot fiber alignment vs distance
            data1 = normal.mean_order_distance;
            data2 = sectio.mean_order_distance;
            confidence1 = normal.mean_order_CI;
            confidence2 = sectio.mean_order_CI;
            data1num = normal.num;
            data2num = sectio.num;
            labelname = 'Local fiber alignment';
            name = 'Histo_order_over_distance';
        case 2
            % plot density vs distance
            data1 = normal.mean_density_distance;
            data2 = sectio.mean_density_distance;
            confidence1 = normal.mean_density_CI;
            confidence2 = sectio.mean_density_CI;
            data1num = normal.num;
            data2num = sectio.num;
            labelname = 'Local fiber area fraction';
            name = 'Histo_density_over_distance'; 
    end

    [fig, nature_data] = make_line_plot_fetal_membranes(bin_edges, data1, data2, confidence1, confidence2, titlename, labelname, paint, alpha, colorname1, colorname2, groupname1, groupname2, [], data1num, data2num);
    
    if j == 1
        ylim([0, 1.1])
    end
    
    % an area fraction above 1 does not exist, same with alignment
    axeshandles = findobj(fig, 'Type', 'axes');
    axeshandles.YTick = axeshandles.YTick(axeshandles.YTick <= 1.05);
   
    exportgraphics(gcf,'', 'Resolution', 600)
    save('', 'nature_data')

    close
end







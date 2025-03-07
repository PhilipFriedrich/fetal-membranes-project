clc, clear, close all


%%
downscale = 2.5;
mum_per_pix =0.25.*downscale;

%%
Annotation1 = dir('');
Annotation2 = dir('');

annotationToUse = Annotation2;

saveTable = [];
saveName = [];
for i = 1:length(annotationToUse)
        tic
        nameToSave = strrep(annotationToUse(i).name, '_mask.tif', '');
        masky = logical(imread([annotationToUse(i).folder, '\', annotationToUse(i).name]));
        masky = imfill(masky, 'holes');
        org = imread(strrep([annotationToUse(i).folder, '\', annotationToUse(i).name], '_mask.tif', '.tif'));
        org_gray = 1- double(rgb2gray(org))./255;
        %rescale
        org = imresize(org, downscale^(-1));
        org_gray = imresize(org_gray, downscale^(-1));
        masky = imresize(masky, downscale^(-1));
        %% color deconvolution
        [Channel_1, Channel_2, Channel_3, stain1, stain2] = ColorDeconvolution_MRE.calcColorTransform(org, 0.99, 1, true, 0.1, 1, false);
        
        hematoxylinSignal = Channel_1 - Channel_2;
        hematoxylinSignal(~masky) = 0;
        hematoxylinMask = imdilate(hematoxylinSignal > 0.1, strel('disk', 3));
        eosinSignal = org_gray; eosinSignal(hematoxylinMask | ~masky | eosinSignal <0) = 0;
        %% calc order
        [S, nx, ny] = get2DQtensor(eosinSignal, mum_per_pix, 5, true,0.0001,0.1, 0.001);

        eosinMask = eosinSignal > otsuthresh(eosinSignal(:))/10;
        toLookAtMask = imerode(masky, strel('disk', round(5*mum_per_pix^(-1)))) & eosinMask;

        S(~toLookAtMask) = NaN;
        %% Thickness of mask
        autoCorr = autocorr2D(masky);
        stepSize = 2;
        maxRange = round(1000*mum_per_pix^(-1));
        [Tics_all, Autocorr2D_all] = calcRadialProfile(autoCorr, stepSize, maxRange);
        Autocorr2D_all(Autocorr2D_all<0) = 0;
        [~, minIndex] = min(abs(exp(-1) - Autocorr2D_all));
        effectiveMaskRadius = Tics_all(minIndex).*mum_per_pix;
        %% Area fractions
        areaFraction_binary = sum(eosinMask, 'all')./sum(masky(:));
        areaFraction_weighted = sum(eosinSignal, 'all')./sum(masky(:));

        averageKernel = fspecial("average", round(5*mum_per_pix^(-1)));
        densityField = conv2(eosinMask, averageKernel, 'same');
        densityField(~toLookAtMask) = NaN;
        %% 
        %only look at biggest object
        CC = bwconncomp(masky);
        [~, maxIndex] = sort(cellfun(@length, CC.PixelIdxList), 'descend');
        masky = zeros(size(masky));
        masky(CC.PixelIdxList{maxIndex(1)}) = 1;
        
        props = regionprops(masky, 'Orientation');
        
        dx = (cos(deg2rad(props.Orientation)));
        dy = (sin(deg2rad(props.Orientation)));

        normalVector = -[-dy, dx]./norm([-dy, dx]);
        translationVector = round(20.*normalVector);

        translatedSegment = false(size(masky));
        
        CC = bwconncomp(masky);
        [xind, yind] = ind2sub(size(masky), CC.PixelIdxList{1});

        for fillIndex = 1:length(xind)
            xTranslated = xind(fillIndex)+translationVector(2);
            yTranslated = yind(fillIndex)-translationVector(1);
            if xTranslated > 0 && xTranslated <= size(masky, 1) && yTranslated > 0 && yTranslated <= size(masky, 2)
                translatedSegment(xTranslated, yTranslated) = true;
            end
        end
        epitheliumBoundary = bwareaopen(((translatedSegment - masky)>0), 15);
        epitheliumBoundary = bwareaopen(epitheliumBoundary, 1000);

        % show epithelium edge with case name
%         figure; imshow(eosinMask);
            

        distanceToEpithelium = bwdist(epitheliumBoundary).*mum_per_pix;
        distances = 10:10:500;
        counter = 0;
        CI_S_distances = [];
        CI_densityField_distances = [];
        distances_NoOfPixels = [];
        correlationFactor = 1 - (round(5*mum_per_pix^(-1))^2-8)/(round(5*mum_per_pix^(-1))^2);
        for distanceIndex = distances
            counter = counter + 1;
            if counter == 1
                currentIndices = distanceToEpithelium < distanceIndex;
                S_distance_mean(counter) = nanmean(S(currentIndices));
                [delta_CI_S] = Calculate_CI_For_UnknownDistribution_correlatedData(S(currentIndices), 0.05, correlationFactor);
                CI_S_distances(counter) = delta_CI_S;
                [delta_CI_density] = Calculate_CI_For_UnknownDistribution_correlatedData(densityField(currentIndices), 0.05, correlationFactor);
                CI_densityField_distances(counter) = delta_CI_density;
                S_distance_std(counter) = nanstd(S(currentIndices));
                density_distance_mean(counter) = nanmean(densityField(currentIndices));
                density_distance_std(counter) = nanstd(densityField(currentIndices));
                distances_NoOfPixels(counter) = correlationFactor.*(sum(currentIndices & ~isnan(densityField), 'all'));
            else
                currentIndices = distanceToEpithelium < distanceIndex & distanceToEpithelium > distances(counter - 1);
                S_distance_mean(counter) = nanmean(S(currentIndices));
                [delta_CI_S] = Calculate_CI_For_UnknownDistribution_correlatedData(S(currentIndices), 0.05, correlationFactor);
                CI_S_distances(counter) = delta_CI_S;
                [delta_CI_density] = Calculate_CI_For_UnknownDistribution_correlatedData(densityField(currentIndices), 0.05, correlationFactor);
                CI_densityField_distances(counter) = delta_CI_density;
                S_distance_std(counter) = nanstd(S(currentIndices));
                density_distance_mean(counter) = nanmean(densityField(currentIndices));
                density_distance_std(counter) = nanstd(densityField(currentIndices));
                distances_NoOfPixels(counter) = correlationFactor.*(sum(currentIndices & ~isnan(densityField), 'all'));
            end
        end
        
%         % show figures
%         figure('Position', [100, 550, 500, 300]); errorbar(distances, density_distance_mean, density_distance_std);  title(nameToSave)   
%         figure('Position', [650, 550, 500, 300]); imagesc(densityField);  title(nameToSave)   
%         figure('Position', [1300, 550, 500, 300]); imshow(epitheliumBoundary); title(nameToSave)   
%         figure('Position', [650, 50, 500, 300]); imshow(org);
%         toc
%         pause()
%         close all

        %save stuff
        saveName{i} = nameToSave;

        saveTable{i, 1} = nanmean(S(:));
        saveTable{i, 2} = nanstd(S(:));
        saveTable{i, 3} = effectiveMaskRadius;
        saveTable{i, 4} = areaFraction_binary;
        saveTable{i, 5} = areaFraction_weighted;
        saveTable{i, 6} = distances;
        saveTable{i, 7} = S_distance_mean;
        saveTable{i, 8} = S_distance_std;
        saveTable{i, 9} = density_distance_mean;
        saveTable{i, 10} = density_distance_std;
        saveTable{i, 11} = distances_NoOfPixels;
        saveTable{i, 12} = CI_densityField_distances;
        saveTable{i, 13} = CI_S_distances;
end

observableNames = {'Mean nematic order parameter', 'STD nematic order parameter',...
    'Length scale of ICT', 'ECM area fraction', 'ECM intensity-weighted area fraction', ...
    'Distances [µm]', 'Mean nematic order distance dependency', 'STD nematic order distance dependency', ...
    'Mean local area fraction (binary, r=5µm) distance dependency', 'STD local area fraction (binary, r=5µm) distance dependency',...
    'Number of pixels considered at each distance', 'Confidence interval density at distance','Confidence interval nematic order at distance' };

if isequal(annotationToUse, Annotation1)
    save('', 'saveTable', 'observableNames', 'saveName');
else
    save('', 'saveTable', 'observableNames', 'saveName');
end
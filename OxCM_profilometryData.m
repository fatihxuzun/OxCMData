function OxCM_profilometryData(nCpu, boundarClean, outlierClean, gx, gy, smth, flipSecondData, outputFilename)
    % STEP 1
    % OxCM_profilometryData processes, aligns, interpolates, and smooths dual-surface
    % profilometry point cloud data (Plane-A and Plane-B) to produce a unified 2.5D 
    % surface profilometry text file for subsequent solid mesh generation.
    %
    % Inputs:
    %   - nCpu: Number of CPU workers for parallel pool (default: 4 or available cores)
    %   - boundarClean: Number of boundary peel iterations to remove perimeter edge points (default: 1)
    %   - outlierClean: Number of outlier removal passes for Z heights (default: 1)
    %   - gx: Common mesh grid spacing in X direction in mm (default: 0.1)
    %   - gy: Common mesh grid spacing in Y direction in mm (default: 0.1)
    %   - smth: B-spline smoothing factor between 0 [fine/less smooth] and 1 [coarse/max smooth] (default: 0.95)
    %   - flipSecondData: Coordinate flipping for plane-B: 'none' (default), 'x', 'y', or 'xy'
    %   - outputFilename: Name of the output text file containing [x, y, z] data (default: 'myData.txt')

    close all; clc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter Defaults & Validation
    if nargin < 1 || isempty(nCpu)
        try
            nCpu = feature('numcores');
        catch
            nCpu = 4;
        end
    end
    if nargin < 2 || isempty(boundarClean)
        boundarClean = 1;
    end
    if nargin < 3 || isempty(outlierClean)
        outlierClean = 1;
    end
    if nargin < 4 || isempty(gx)
        gx = 0.1;
    end
    if nargin < 5 || isempty(gy)
        gy = gx;
    end
    if nargin < 6 || isempty(smth)
        smth = 0.95;
    end
    if nargin < 7 || isempty(flipSecondData)
        flipSecondData = 'none';
    end
    if nargin < 8 || isempty(outputFilename)
        outputFilename = 'myData.txt';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if a parallel pool is currently active
    parCheck = gcp('nocreate');
    
    % If a pool is active but has a different size, delete it
    if ~isempty(parCheck) && parCheck.NumWorkers ~= nCpu
        delete(parCheck);
        parCheck = [];
    end
    
    % If no pool is active, start a new one with 'nCpu' workers
    if isempty(parCheck)
        myCluster = parcluster('local');
        myCluster.NumWorkers = nCpu;
        parpool(myCluster, nCpu);
    end
    disp('Parallel pool is ACTIVE');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read data files
    disp('Select plane-A profilometry data file:');
    [fileI, pathI] = uigetfile('*.txt', 'Select plane-A profilometry data file');
    if isequal(fileI, 0)
        warning('File selection cancelled for Plane-A. Aborting.');
        return;
    end

    disp('Select plane-B profilometry data file:');
    [fileJ, pathJ] = uigetfile('*.txt', 'Select plane-B profilometry data file');
    if isequal(fileJ, 0)
        warning('File selection cancelled for Plane-B. Aborting.');
        return;
    end
    
    rawA = importdata(fullfile(pathI, fileI), ' ');
    if isstruct(rawA), dataA = rawA.data'; else, dataA = rawA'; end

    rawB = importdata(fullfile(pathJ, fileJ), ' ');
    if isstruct(rawB), dataB = rawB.data'; else, dataB = rawB'; end
    
    % Move data centers to the origin
    dataA(3,:) = dataA(3,:) - mean(dataA(3,:));
    dataB(3,:) = dataB(3,:) - mean(dataB(3,:));
    dataA(2,:) = dataA(2,:) - mean(dataA(2,:));
    dataB(2,:) = dataB(2,:) - mean(dataB(2,:));
    dataA(1,:) = dataA(1,:) - mean(dataA(1,:));
    dataB(1,:) = dataB(1,:) - mean(dataB(1,:));

    % Flip plane-B coordinates if requested
    switch lower(flipSecondData)
        case 'x'
            dataB(1,:) = -dataB(1,:);
        case 'y'
            dataB(2,:) = -dataB(2,:);
        case 'xy'
            dataB(1:2,:) = -dataB(1:2,:);
        case 'none'
            % No flipping applied
        otherwise
            warning('Invalid flipSecondData option ("%s"). No flipping applied.', flipSecondData);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean the data at out-of-bounds (boundaries)
    for l = 1:boundarClean
        tic
        kA = boundary(dataA(1,:)', dataA(2,:)');
        dataA(:,kA) = [];  
        
        kB = boundary(dataB(1,:)', dataB(2,:)');
        dataB(:,kB) = [];
        toc
    end

    % Clean outlier points along the Z axis
    for l = 1:outlierClean
        [~, outA] = rmoutliers(dataA(3,:));
        dataA(:, outA == 1) = [];
        [~, outB] = rmoutliers(dataB(3,:));
        dataB(:, outB == 1) = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create XY reference plane data for flattening
    icpDataLen = 10000;
    
    xmax = max(dataA(1,:)) / 2;
    xmin = min(dataA(1,:)) / 2;
    ymax = max(dataA(2,:)) / 2;
    ymin = min(dataA(2,:)) / 2;
    
    algxy = zeros(3, icpDataLen);
    algxy(1,:) = xmin + (xmax - xmin) * rand(1, icpDataLen);
    algxy(2,:) = ymin + (ymax - ymin) * rand(1, icpDataLen);
    
    rSize = round(length(dataA) / icpDataLen);
    if rSize < 1
        rSize = 1;
    end
    
    algA = dataA(1:3, 1:rSize:end);
    algB = dataB(1:3, 1:rSize:end);
    
    [TR, TT] = hDICIcp(algxy, algA, 10);
    algA = TR * algA + repmat(TT, 1, length(algA));
    dataA = TR * dataA + repmat(TT, 1, length(dataA));
    
    [TR, TT] = hDICIcp(algxy, algB, 10);
    algB = TR * algB + repmat(TT, 1, length(algB));
    dataB = TR * dataB + repmat(TT, 1, length(dataB));

    [TR, TT] = hDICIcp(algA, algB, 50);
    dataB = TR * dataB + repmat(TT, 1, length(dataB));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-center data in XY plane
    dataA(2,:) = dataA(2,:) - mean(dataA(2,:));
    dataB(2,:) = dataB(2,:) - mean(dataB(2,:));
    dataA(1,:) = dataA(1,:) - mean(dataA(1,:));
    dataB(1,:) = dataB(1,:) - mean(dataB(1,:));

    d = mean(abs(dataA(3,:))) * 5 / ((mean(abs(dataA(1,:))) + mean(abs(dataA(2,:)))) / 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1); 
    scatter3(dataA(1,:), dataA(2,:), dataA(3,:), 1, '.', 'r'); hold on; 
    scatter3(dataB(1,:), dataB(2,:), dataB(3,:), 1, '.', 'b'); hold off; 
    shading interp; view(60,15); grid off; axis equal; axis tight; daspect([1, 1, d]); title('Aligned');
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); legend('First','Second','Location','east'); 
    fontsize(16, 'points'); fontname('Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xMinR = min(dataA(1,:)); xMaxR = max(dataA(1,:));
    yMinR = min(dataA(2,:)); yMaxR = max(dataA(2,:));
    coordH = [xMinR, yMaxR; xMaxR, yMaxR; xMaxR, yMinR; xMinR, yMinR];

    % Determine meshgrid density
    xmax = max(coordH(:,1)) + gx;
    xmin = min(coordH(:,1));
    ymax = max(coordH(:,2)) + gy;
    ymin = min(coordH(:,2));
    [x, y] = meshgrid(xmin:gx:xmax, ymax:-gy:ymin);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interpolate to common meshgrid
    kA = boundary(dataA(1,:)', dataA(2,:)');
    inA = inpolygon(x, y, dataA(1,kA), dataA(2,kA));
    zi = griddata(dataA(1,:), dataA(2,:), dataA(3,:), x, y, 'cubic'); 
    zi(~inA) = NaN;

    kB = boundary(dataB(1,:)', dataB(2,:)');
    inB = inpolygon(x, y, dataB(1,kB), dataB(2,kB));
    zj = griddata(dataB(1,:), dataB(2,:), dataB(3,:), x, y, 'cubic'); 
    zj(~inB) = NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot interpolated figures
    figure(2); 
    surf(x, y, zi); hold on; surf(x, y, zj); hold off; 
    shading interp; view(60,15); grid off; colormap jet; axis equal; axis tight; daspect([1, 1, d]); 
    title('Interpolated to Common Grid');
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); fontsize(16, 'points'); fontname('Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Average and spline-smooth the surface
    zTemp = (zi + zj) / 2;
    zFilled = fillmissing2(zTemp, 'nearest');

    c = {y(:, 1)', x(1, :)};

    numYPieces = max(1, round(1 + (size(zTemp, 1) - 4) * (1 - smth)));
    numXPieces = max(1, round(1 + (size(zTemp, 2) - 4) * (1 - smth)));

    smthfunc = spap2({numYPieces, numXPieces}, [4 4], c, zFilled);

    zSmth = fnval(smthfunc, c);
    zSmth(isnan(zTemp)) = NaN;

    figure(3); 
    surf(x, y, zTemp); shading interp; view(60,15); grid off; colormap jet; axis equal; axis tight; daspect([1, 1, d]); 
    title('Averaged Raw'); xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); fontsize(16, 'points'); fontname('Arial');

    figure(4); 
    surf(x, y, zSmth); shading interp; view(60,15); grid off; colormap jet; axis equal; axis tight; daspect([1, 1, d]); 
    title('Smoothed'); xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis'); fontsize(16, 'points'); fontname('Arial');

    figure(5); 
    plot(mean(x, 'omitnan'), mean(zTemp, 'omitnan'), '-r'); hold on;
    plot(mean(x, 'omitnan'), mean(zSmth, 'omitnan'), '-b'); 
    title('Averaged (Red) vs Smoothed (Blue) Data'); 
    xlabel('x-axis'); ylabel('z-axis'); legend('Averaged Raw','Smoothed','Location','east'); 
    fontsize(16, 'points'); fontname('Arial');

    figure(6); 
    plot(mean(y', 'omitnan'), mean(zTemp, 'omitnan'), '-r'); hold on;
    plot(mean(y', 'omitnan'), mean(zSmth, 'omitnan'), '-b'); 
    title('Averaged (Red) vs Smoothed (Blue) Data'); 
    xlabel('y-axis'); ylabel('z-axis'); legend('Averaged Raw','Smoothed','Location','east'); 
    fontsize(16, 'points'); fontname('Arial');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write output point cloud to file
    validIdx = ~isnan(zTemp);
    myData = [x(validIdx), y(validIdx), zSmth(validIdx)];
    myData = single(myData);
    
    writematrix(myData, outputFilename, 'Delimiter', '\t');
    fprintf('Successfully written output profilometry data to "%s".\n', outputFilename);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Custom Parallel Iterative Closest Point (ICP) Algorithm
    function [rFinal, tFinal, rMSErrors, runTimes] = hDICIcp(refCloud, movCloud, varargin)
        p = inputParser;
        p.addRequired('refCloud', @(x) isreal(x) && size(x,1) == 3);
        p.addRequired('movCloud', @(x) isreal(x) && size(x,1) == 3);
        p.addOptional('maxters', 10, @(x) x > 0 && x < 10^5);
        p.addParameter('BndIdx', [], @(x) size(x,1) == 1);
        p.addParameter('FilterEdges', false, @(x) islogical(x));
        p.addParameter('UseExtrap', false, @(x) islogical(x));
        p.addParameter('SearchType', 'bruteForce', @(x) any(strcmpi(x, {'bruteForce','Delaunay','kDtree'})));
        p.addParameter('MinType', 'point', @(x) any(strcmpi(x, {'point','plane','lmapoint'})));
        p.addParameter('NormalsIn', [], @(x) isreal(x) && size(x,1) == 3);
        p.addParameter('SaveAllSteps', false, @(x) islogical(x));
        p.addParameter('MeshTri', [], @(x) isreal(x) && size(x,2) == 3);
        p.addParameter('Logging', false, @(x) islogical(x));
        p.addParameter('WeightFunc', @(x) ones(1,length(x)), @(x) isa(x,'function_handle'));
        p.addParameter('OutlierPct', 0, @(x) isscalar(x) && x > 0 && x < 1);
        
        p.parse(refCloud, movCloud, varargin{:});
        opts = p.Results;
        
        runTimes = zeros(opts.maxters + 1, 1); 
        tic;
        
        numMovPts = size(movCloud, 2);
        currCloud = movCloud;
        rMSErrors = zeros(opts.maxters + 1, 1); 
        
        tTmp = zeros(3, 1);
        rTmp = eye(3);
        
        tHistory = zeros(3, 1, opts.maxters + 1);
        rHistory = repmat(eye(3), [1, 1, opts.maxters + 1]);
            
        if strcmp(opts.MinType, 'plane') && isempty(opts.NormalsIn)
            opts.NormalsIn = estimateNormalsParallel(refCloud, 4);
        end
        
        if strcmp(opts.SearchType, 'Delaunay')
            delaunayObj = DelaunayTri(transpose(refCloud));
        end
        
        if strcmp(opts.SearchType, 'kDtree')
            kdTreeObj = KDTreeSearcher(transpose(refCloud));
        end
        
        if opts.FilterEdges
            if isempty(opts.BndIdx)
                edgeNodes = extractBoundary(refCloud, opts.MeshTri);
            else
                edgeNodes = opts.BndIdx;
            end
        end
        
        if opts.UseExtrap
            quatHist = [ones(1, opts.maxters + 1); zeros(6, opts.maxters + 1)];   
            deltaQuat = zeros(7, opts.maxters + 1);
            angChange = zeros(1, opts.maxters + 1);
        end
        
        runTimes(1) = toc;
        
        for iterIdx = 1:opts.maxters
            switch opts.SearchType
                case 'bruteForce'
                    [pairs, distances] = matchBruteForceParallel(refCloud, currCloud);
                case 'Delaunay'
                    [pairs, distances] = matchDelaunay(refCloud, currCloud, delaunayObj);
                case 'kDtree'
                    [pairs, distances] = matchKdTree(currCloud, kdTreeObj);
            end
        
            if opts.FilterEdges
                validMask = ~ismember(pairs, edgeNodes);
                refIdx = pairs(validMask);
                distances = distances(validMask);
            else
                validMask = true(1, numMovPts);
                refIdx = pairs;
            end
            
            if opts.OutlierPct
                cutoff = round((1 - opts.OutlierPct) * sum(validMask));
                activePairs = find(validMask);
                [~, sortedIdx] = sort(distances);
                validMask(activePairs(sortedIdx(cutoff:end))) = false;
                refIdx = pairs(validMask);
                distances = distances(validMask);
            end
            
            if iterIdx == 1
                rMSErrors(iterIdx) = sqrt(sum(distances.^2) / length(distances));
            end
            
            switch lower(opts.MinType)
                case 'point'
                    wVector = opts.WeightFunc(pairs);
                    [rTmp, tTmp] = solvePointMin(refCloud(:, refIdx), currCloud(:, validMask), wVector(validMask));
                case 'plane'
                    wVector = opts.WeightFunc(pairs);
                    [rTmp, tTmp] = solvePlaneMin(refCloud(:, refIdx), currCloud(:, validMask), opts.NormalsIn(:, refIdx), wVector(validMask));
                case 'lmapoint'
                    [rTmp, tTmp] = solveLmaMin(refCloud(:, refIdx), currCloud(:, validMask));
            end
        
            rHistory(:, :, iterIdx + 1) = rTmp * rHistory(:, :, iterIdx);
            tHistory(:, :, iterIdx + 1) = rTmp * tHistory(:, :, iterIdx) + tTmp;
            
            currCloud = rHistory(:, :, iterIdx + 1) * movCloud + repmat(tHistory(:, :, iterIdx + 1), 1, numMovPts);
            rMSErrors(iterIdx + 1) = calculateRms(refCloud(:, refIdx), currCloud(:, validMask));
            
            if opts.UseExtrap
                quatHist(:, iterIdx + 1) = [matrixToQuat(rHistory(:, :, iterIdx + 1)); tHistory(:, :, iterIdx + 1)];
                deltaQuat(:, iterIdx + 1) = quatHist(:, iterIdx + 1) - quatHist(:, iterIdx);
                angChange(iterIdx + 1) = (180 / pi) * acos(dot(deltaQuat(:, iterIdx), deltaQuat(:, iterIdx + 1)) / (norm(deltaQuat(:, iterIdx)) * norm(deltaQuat(:, iterIdx + 1))));
                
                if opts.Logging
                    disp(['Direction change: ' num2str(angChange(iterIdx + 1)) ' deg. Iteration: ' num2str(iterIdx)]);
                end
                
                if iterIdx > 2 && angChange(iterIdx + 1) < 10 && angChange(iterIdx) < 10
                    eVal = [rMSErrors(iterIdx + 1), rMSErrors(iterIdx), rMSErrors(iterIdx - 1)];
                    vVal = [0, -norm(deltaQuat(:, iterIdx + 1)), -norm(deltaQuat(:, iterIdx)) - norm(deltaQuat(:, iterIdx + 1))];
                    vLimit = 25 * norm(deltaQuat(:, iterIdx + 1));
                    stepSize = performExtrapolation(vVal, eVal, vLimit);
                    
                    if stepSize ~= 0
                        predictedQuat = quatHist(:, iterIdx + 1) + stepSize * deltaQuat(:, iterIdx + 1) / norm(deltaQuat(:, iterIdx + 1));
                        predictedQuat(1:4) = predictedQuat(1:4) / norm(predictedQuat(1:4));
                        quatHist(:, iterIdx + 1) = predictedQuat;
                        rHistory(:, :, iterIdx + 1) = quatToMatrix(quatHist(1:4, iterIdx + 1));
                        tHistory(:, :, iterIdx + 1) = quatHist(5:7, iterIdx + 1);
                        
                        currCloud = rHistory(:, :, iterIdx + 1) * movCloud + repmat(tHistory(:, :, iterIdx + 1), 1, numMovPts);
                        
                        switch opts.SearchType
                            case 'bruteForce'
                                [~, distances] = matchBruteForceParallel(refCloud, currCloud);
                            case 'Delaunay'
                                [~, distances] = matchDelaunay(refCloud, currCloud, delaunayObj);
                            case 'kDtree'
                                [~, distances] = matchKdTree(currCloud, kdTreeObj);
                        end
                        rMSErrors(iterIdx + 1) = sqrt(sum(distances.^2) / length(distances));
                    end
                end
            end
            runTimes(iterIdx + 1) = toc;
        end
        
        if ~opts.SaveAllSteps
            rFinal = rHistory(:, :, end);
            tFinal = tHistory(:, :, end);
        else
            rFinal = rHistory;
            tFinal = tHistory;
        end
    end
    
    %% --- Parallelized Subfunctions ---
    function [matchIdx, minDists] = matchBruteForceParallel(ref, mov)
        m = size(mov, 2);
        n = size(ref, 2);    
        matchIdx = zeros(1, m);
        minDists = zeros(1, m);
        
        parfor i = 1:m
            sqDists = zeros(1, n);
            ptMov = mov(:, i);
            for dim = 1:3
                sqDists = sqDists + (ref(dim, :) - ptMov(dim)).^2;
            end
            [minVal, minLoc] = min(sqDists);
            minDists(i) = minVal;
            matchIdx(i) = minLoc;
        end
        minDists = sqrt(minDists);
    end
    
    function nVecs = estimateNormalsParallel(pts, kNeighbors)
        totalPts = size(pts, 2);
        nVecs = zeros(3, totalPts);
        toolboxCheck = ver('stats');
        
        if ~isempty(toolboxCheck) && str2double(toolboxCheck.Version) >= 7.5 
            nbIds = transpose(knnsearch(transpose(pts), transpose(pts), 'k', kNeighbors + 1));
        else
            nbIds = fallbackKnnParallel(pts, pts, kNeighbors + 1);
        end
        
        parfor i = 1:totalPts
            neighborhood = pts(:, nbIds(2:end, i));
            localMean = mean(neighborhood, 2);
            centered = neighborhood - localMean;
            covMat = centered * centered';
            [V, D] = eig(covMat);
            [~, minIdx] = min(diag(D));
            nVecs(:, i) = V(:, minIdx);   
        end
    end
    
    function [nbIds, nbDists] = fallbackKnnParallel(dataMat, queryMat, kNeighbors)
        numQueries = size(queryMat, 2);
        numData = size(dataMat, 2);
        nbIds = zeros(kNeighbors, numQueries);
        nbDists = zeros(kNeighbors, numQueries);
        dims = size(dataMat, 1);
        
        parfor i = 1:numQueries
            localDists = zeros(1, numData);
            qPt = queryMat(:, i);
            for dim = 1:dims
                localDists = localDists + (dataMat(dim, :) - qPt(dim)).^2;
            end
            
            tempIds = zeros(kNeighbors, 1);
            tempDists = zeros(kNeighbors, 1);
            for j = 1:kNeighbors
                [val, idx] = min(localDists);
                tempIds(j) = idx;
                tempDists(j) = sqrt(val);
                localDists(idx) = NaN;
            end
            nbIds(:, i) = tempIds;
            nbDists(:, i) = tempDists;
        end
    end
    
    %% --- Vectorized/Standard Subfunctions ---
    function [matchIdx, minDists] = matchDelaunay(ref, mov, delaunayObj)
        matchIdx = transpose(nearestNeighbor(delaunayObj, transpose(mov)));
        minDists = sqrt(sum((mov - ref(:, matchIdx)).^2, 1));
    end
    
    function [matchIdx, minDists] = matchKdTree(mov, kdTreeObj)
        [matchIdx, minDists] = knnsearch(kdTreeObj, transpose(mov));
        matchIdx = transpose(matchIdx);
    end
    
    function [R, T] = solvePointMin(ref, mov, weights)
        wNorm = weights ./ sum(weights);
        refCentroid = ref * wNorm';
        refDev = (ref - refCentroid) .* wNorm;
        movCentroid = mov * wNorm';
        movDev = mov - movCentroid;
        
        S = movDev * refDev'; 
        [U, ~, V] = svd(S);
        R = V * diag([1 1 det(U * V')]) * U';
        T = refCentroid - R * movCentroid;
    end
    
    function [R, T] = solvePlaneMin(ref, mov, normals, weights)
        scaledNormals = normals .* weights;
        crossProducts = cross(mov, scaledNormals);
        combinedGeom = vertcat(crossProducts, scaledNormals);
        covGeom = combinedGeom * combinedGeom';
        
        diffVec = mov - ref;
        residualB = - [sum(sum(diffVec .* combinedGeom(1,:) .* normals));
                       sum(sum(diffVec .* combinedGeom(2,:) .* normals));
                       sum(sum(diffVec .* combinedGeom(3,:) .* normals));
                       sum(sum(diffVec .* combinedGeom(4,:) .* normals));
                       sum(sum(diffVec .* combinedGeom(5,:) .* normals));
                       sum(sum(diffVec .* combinedGeom(6,:) .* normals))];
           
        paramX = covGeom \ residualB;
        cx = cos(paramX(1)); cy = cos(paramX(2)); cz = cos(paramX(3)); 
        sx = sin(paramX(1)); sy = sin(paramX(2)); sz = sin(paramX(3)); 
        
        R = [cy*cz cz*sx*sy-cx*sz cx*cz*sy+sx*sz;
             cy*sz cx*cz+sx*sy*sz cx*sy*sz-cz*sx;
             -sy cy*sx cx*cy];
        T = paramX(4:6);
    end
    
    function [R, T] = solveLmaMin(ref, mov)
        rx = @(a) [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
        ry = @(b) [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
        rz = @(g) [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];
        getRot = @(x) rx(x(1)) * ry(x(2)) * rz(x(3));
        mapFunc = @(x, data) getRot(x(1:3)) * data + repmat(x(4:6), 1, length(data));
        
        cfg = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off');
        solvedVals = lsqcurvefit(mapFunc, zeros(6,1), mov, ref, [], [], cfg);
        
        R = getRot(solvedVals(1:3));
        T = solvedVals(4:6);
    end
    
    function step = performExtrapolation(v, d, maxistep)
        p1 = polyfit(v, d, 1); 
        p2 = polyfit(v, d, 2); 
        v1 = -p1(2) / p1(1); 
        v2 = -p2(2) / (2 * p2(1)); 
        
        if issorted([0 v2 v1 maxistep]) || issorted([0 v2 maxistep v1])
            step = v2;
        elseif issorted([0 v1 v2 maxistep]) || issorted([0 v1 maxistep v2]) || (v2 < 0 && issorted([0 v1 maxistep]))
            step = v1;
        elseif v1 > maxistep && v2 > maxistep
            step = maxistep;
        else
            step = 0;
        end
    end
    
    function err = calculateRms(p1, p2)
        err = sqrt(mean(sum((p1 - p2).^2, 1)));
    end
    
    function quat = matrixToQuat(R)
        w = 0.5 * sqrt(1 + R(1,1) + R(2,2) + R(3,3));
        x = 0.5 * sign(R(3,2) - R(2,3)) * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        y = 0.5 * sign(R(1,3) - R(3,1)) * sqrt(1 - R(1,1) + R(2,2) - R(3,3));
        z = 0.5 * sign(R(2,1) - R(1,2)) * sqrt(1 - R(1,1) - R(2,2) + R(3,3));
        quat = reshape([w; x; y; z], 4, []);
    end
    
    function R = quatToMatrix(q)
        q0(1,1,:) = q(1,:); qx(1,1,:) = q(2,:); qy(1,1,:) = q(3,:); qz(1,1,:) = q(4,:);
        R = [q0.^2+qx.^2-qy.^2-qz.^2 2*qx.*qy-2*q0.*qz 2*qx.*qz+2*q0.*qy;
             2*qx.*qy+2*q0.*qz q0.^2-qx.^2+qy.^2-qz.^2 2*qy.*qz-2*q0.*qx;
             2*qx.*qz-2*q0.*qy 2*qy.*qz+2*q0.*qx q0.^2-qx.^2-qy.^2+qz.^2];
    end
    
    function boundaryNodes = extractBoundary(pts, triMatrix)
        TR = TriRep(double(triMatrix), double(pts(1,:))', double(pts(2,:))', double(pts(3,:))');
        FF = freeBoundary(TR);
        boundaryNodes = FF(:, 1);
    end
end
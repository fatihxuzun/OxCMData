function OxCM_extrudedPerimeterMesh(filename, extrusionDepth, hmax, hmin, scaleFactor, numLayersZ, outputFilename, elementType)
    % STEP 3
    % OxCM_extrudedPerimeterMesh generates an extruded 3D solid mesh composed of 
    % parallel structured Z-layers from 2.5D profilometry point cloud data.
    % 
    % Inputs:
    % - filename: Name of input text file containing [x, y, z] data (default: 'myData.txt')
    % - extrusionDepth: Total Z extrusion depth in mm (default: 10)
    % - hmax: Target maximum 2D element edge length in XY plane (default: [] -> raw cloud)
    % - hmin: Target minimum 2D element edge length in XY plane (default: [] -> raw cloud)
    % - scaleFactor: Factor to shrink the 2D plane inwards (default: 1.0)
    % - numLayersZ: Explicit number of parallel Z-layers (default: 10)
    % - outputFilename: Name of the output XML file (default: 'myMesh.xml')
    % - elementType: 'tetrahedron' (default) or 'hexahedron'
    
    % --- Input Handling ---
    if nargin < 1 || isempty(filename), filename = 'myData.txt'; end
    if nargin < 2 || isempty(extrusionDepth), extrusionDepth = 10; end
    if nargin < 3, hmax = []; end
    if nargin < 4, hmin = []; end 
    if nargin < 5 || isempty(scaleFactor), scaleFactor = 1.0; end 
    if nargin < 6 || isempty(numLayersZ), numLayersZ = 10; end 
    if nargin < 7 || isempty(outputFilename), outputFilename = 'myMesh.xml'; end
    if nargin < 8 || isempty(elementType), elementType = 'tetrahedron'; end
    
    isHex = strncmpi(elementType, 'hex', 3);
    
    if ~isfile(filename)
        error('Input file "%s" not found. Please check the path.', filename);
    end
    data = readmatrix(filename);
    x = data(:,1); y = data(:,2);
    
    % --- Centroid-Based 2D Scaling ---
    if scaleFactor ~= 1.0
        cx = mean(x); cy = mean(y);
        x = (x - cx) * scaleFactor + cx;
        y = (y - cy) * scaleFactor + cy;
    end
    
    shp = alphaShape(x, y);
    shp.Alpha = shp.Alpha * 1.2; 
    top_tris = alphaTriangulation(shp);
    
    % Clean unreferenced nodes
    used_nodes = unique(top_tris(:));
    node_map = zeros(max(used_nodes), 1);
    node_map(used_nodes) = 1:length(used_nodes);
    x = x(used_nodes); y = y(used_nodes);
    top_tris = node_map(top_tris);
    
    % --- 2D Surface Remeshing with hmax / hmin ---
    if (~isempty(hmax) || ~isempty(hmin)) && ~isempty(ver('pde'))
        try
            model2d = createpde();
            geometryFromMesh(model2d, [x'; y'], top_tris');
            
            meshArgs = {'GeometricOrder', 'linear'};
            if ~isempty(hmax), meshArgs = [meshArgs, {'Hmax', hmax}]; end
            if ~isempty(hmin), meshArgs = [meshArgs, {'Hmin', hmin}]; end
            
            generateMesh(model2d, meshArgs{:});
            x = model2d.Mesh.Nodes(1,:)';
            y = model2d.Mesh.Nodes(2,:)';
            top_tris = model2d.Mesh.Elements';
        catch ME
            warning('2D Remeshing failed (%s). Falling back to raw point cloud triangulation.', ME.message);
        end
    end
    
    num_layers = max(1, round(numLayersZ));
    z_vals = linspace(0, -abs(extrusionDepth), num_layers + 1);

    %% ======================= 1. HEXAHEDRAL BRANCH ======================= %%
    if isHex
        % Ensure 2D triangles are oriented Counter-Clockwise (CCW)
        v12 = [x(top_tris(:,2)) - x(top_tris(:,1)), y(top_tris(:,2)) - y(top_tris(:,1))];
        v13 = [x(top_tris(:,3)) - x(top_tris(:,1)), y(top_tris(:,3)) - y(top_tris(:,1))];
        area2 = v12(:,1).*v13(:,2) - v12(:,2).*v13(:,1);
        flip_idx = area2 < 0;
        top_tris(flip_idx, [2 3]) = top_tris(flip_idx, [3 2]);
        
        % Tri-to-Quad Barycentric Subdivision:
        % Find unique edges and their midpoints
        Nt = size(top_tris, 1);
        Nv = length(x);
        all_edges = [top_tris(:, [1 2]); top_tris(:, [2 3]); top_tris(:, [3 1])];
        sorted_edges = sort(all_edges, 2);
        [u_edges, ~, edge_map] = unique(sorted_edges, 'rows');
        Ne = size(u_edges, 1);
        
        edge_mid_x = (x(u_edges(:,1)) + x(u_edges(:,2))) / 2;
        edge_mid_y = (y(u_edges(:,1)) + y(u_edges(:,2))) / 2;
        
        % Centroids of triangles
        cent_x = (x(top_tris(:,1)) + x(top_tris(:,2)) + x(top_tris(:,3))) / 3;
        cent_y = (y(top_tris(:,1)) + y(top_tris(:,2)) + y(top_tris(:,3))) / 3;
        
        % Combined 2D Nodes
        x_2d = [x; edge_mid_x; cent_x];
        y_2d = [y; edge_mid_y; cent_y];
        N = length(x_2d); % Total 2D nodes per Z-layer
        
        % Split each triangle into 3 CCW Quadrilaterals
        v1 = top_tris(:, 1);
        v2 = top_tris(:, 2);
        v3 = top_tris(:, 3);
        m1 = Nv + edge_map(1:Nt);
        m2 = Nv + edge_map(Nt+1:2*Nt);
        m3 = Nv + edge_map(2*Nt+1:3*Nt);
        c  = Nv + Ne + (1:Nt)';
        
        top_quads = [v1, m1, c, m3;
                     v2, m2, c, m1;
                     v3, m3, c, m2];
                 
        % Build 3D Nodes Array
        Nodes = zeros(N * (num_layers + 1), 3);
        for k = 1:(num_layers + 1)
            idx_start = (k - 1) * N + 1;
            idx_end = k * N;
            Nodes(idx_start:idx_end, :) = [x_2d, y_2d, repmat(z_vals(k), N, 1)];
        end
        
        % Extrude 2D Quads into 3D Hexahedra
        num_quads = size(top_quads, 1);
        total_hex = num_quads * num_layers;
        Elements = zeros(total_hex, 8);
        
        hex_idx = 1;
        for k = 1:num_layers
            % Bottom face at lower Z (layer k+1), Top face at higher Z (layer k)
            O_bot = k * N;       
            O_top = (k - 1) * N; 
            
            b1 = top_quads(:,1) + O_bot;
            b2 = top_quads(:,2) + O_bot;
            b3 = top_quads(:,3) + O_bot;
            b4 = top_quads(:,4) + O_bot;
            
            t1 = top_quads(:,1) + O_top;
            t2 = top_quads(:,2) + O_top;
            t3 = top_quads(:,3) + O_top;
            t4 = top_quads(:,4) + O_top;
            
            % FEniCS Hex Canonical Ordering: [b1, b2, b3, b4, t1, t2, t3, t4]
            Elements(hex_idx : hex_idx + num_quads - 1, :) = [b1, b2, b3, b4, t1, t2, t3, t4];
            hex_idx = hex_idx + num_quads;
        end
        
        % Ensure Positive Jacobian Determinant for Hexahedra
        P1 = Nodes(Elements(:,1), :);
        P2 = Nodes(Elements(:,2), :);
        P4 = Nodes(Elements(:,4), :);
        P5 = Nodes(Elements(:,5), :);
        
        V12 = P2 - P1; V14 = P4 - P1; V15 = P5 - P1;
        crossV = [V12(:,2).*V14(:,3) - V12(:,3).*V14(:,2), ...
                  V12(:,3).*V14(:,1) - V12(:,1).*V14(:,3), ...
                  V12(:,1).*V14(:,2) - V12(:,2).*V14(:,1)];
        detV = sum(crossV .* V15, 2);
        
        negIdx = detV < 0;
        if any(negIdx)
            Elements(negIdx, [2 4 6 8]) = Elements(negIdx, [4 2 8 6]);
        end
        
        cellType = 'hexahedron';

    %% ===================== 2. TETRAHEDRAL BRANCH ====================== %%
    else
        N = length(x);
        Nodes = zeros(N * (num_layers + 1), 3);
        for k = 1:(num_layers + 1)
            idx_start = (k - 1) * N + 1;
            idx_end = k * N;
            Nodes(idx_start:idx_end, :) = [x, y, repmat(z_vals(k), N, 1)];
        end
        
        sorted_tris = sort(top_tris, 2);
        num_tris = size(sorted_tris, 1);
        num_tets_per_layer = num_tris * 3;
        total_tets = num_tets_per_layer * num_layers;
        Elements = zeros(total_tets, 4);
        
        tet_idx = 1;
        for k = 1:num_layers
            O0 = (k - 1) * N;
            O1 = k * N;
            
            u1 = sorted_tris(:,1) + O0;
            u2 = sorted_tris(:,2) + O0;
            u3 = sorted_tris(:,3) + O0;
            
            v1 = sorted_tris(:,1) + O1;
            v2 = sorted_tris(:,2) + O1;
            v3 = sorted_tris(:,3) + O1;
            
            tets_k = [
                u1, u2, u3, v3;
                u1, u2, v3, v2;
                u1, v2, v3, v1
            ];
            
            Elements(tet_idx : tet_idx + num_tets_per_layer - 1, :) = tets_k;
            tet_idx = tet_idx + num_tets_per_layer;
        end
        
        % Ensure Positive Jacobian Determinant for Tetrahedra
        P1 = Nodes(Elements(:,1), :);
        P2 = Nodes(Elements(:,2), :);
        P3 = Nodes(Elements(:,3), :);
        P4 = Nodes(Elements(:,4), :);
        
        V12 = P2 - P1; V13 = P3 - P1; V14 = P4 - P1;
        crossV = [V12(:,2).*V13(:,3) - V12(:,3).*V13(:,2), ...
                  V12(:,3).*V14(:,1) - V12(:,1).*V13(:,3), ...
                  V12(:,1).*V13(:,2) - V12(:,2).*V13(:,1)];
        detV = sum(crossV .* V14, 2);
        
        negIdx = detV < 0;
        if any(negIdx)
            Elements(negIdx, [1 2]) = Elements(negIdx, [2 1]);
        end
        
        cellType = 'tetrahedron';
    end
    
    % --- Export to FEniCS XML ---
    writeFenicsXML(Nodes, Elements, outputFilename, cellType);
    disp(['Total Nodes: ', num2str(size(Nodes, 1)), ' | Total 3D Elements: ', num2str(size(Elements, 1)), ' (', cellType, ')']);
end

%% ========================== HELPER EXPORTER ========================== %%
function writeFenicsXML(pts, conn, xmlFilename, cellType)
    fid = fopen(xmlFilename, 'w');
    if fid == -1
        error('Cannot open or create file "%s". Check write permissions or folder path.', xmlFilename);
    end
    
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<dolfin xmlns:dolfin="http://fenicsproject.org">\n');
    fprintf(fid, '  <mesh celltype="%s" dim="3">\n', cellType);
    
    % Export Vertices
    fprintf(fid, '    <vertices size="%d">\n', size(pts, 1));
    vData = [0:(size(pts,1)-1); pts(:,1)'; pts(:,2)'; pts(:,3)'];
    fprintf(fid, '      <vertex index="%d" x="%f" y="%f" z="%f" />\n', vData);
    fprintf(fid, '    </vertices>\n');
    
    % Export Cells
    fprintf(fid, '    <cells size="%d">\n', size(conn, 1));
    if strcmp(cellType, 'tetrahedron')
        cData = [0:(size(conn,1)-1); conn(:,1)'-1; conn(:,2)'-1; conn(:,3)'-1; conn(:,4)'-1];
        fprintf(fid, '      <tetrahedron index="%d" v0="%d" v1="%d" v2="%d" v3="%d" />\n', cData);
    elseif strcmp(cellType, 'hexahedron')
        cData = [0:(size(conn,1)-1); conn(:,1)'-1; conn(:,2)'-1; conn(:,3)'-1; conn(:,4)'-1; ...
                                     conn(:,5)'-1; conn(:,6)'-1; conn(:,7)'-1; conn(:,8)'-1];
        fprintf(fid, '      <hexahedron index="%d" v0="%d" v1="%d" v2="%d" v3="%d" v4="%d" v5="%d" v6="%d" v7="%d" />\n', cData);
    else
        cData = [0:(size(conn,1)-1); conn(:,1)'-1; conn(:,2)'-1; conn(:,3)'-1];
        fprintf(fid, '      <triangle index="%d" v0="%d" v1="%d" v2="%d" />\n', cData);
    end
    fprintf(fid, '    </cells>\n  </mesh>\n</dolfin>\n');
    fclose(fid);
end
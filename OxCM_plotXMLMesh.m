function OxCM_plotXMLMesh(xmlFilename)
    % STEP 4
    % OxCM_plotXMLMesh reads a FEniCS/DOLFIN XML mesh file, automatically
    % detects the cell type (tetrahedron, hexahedron, triangle), and plots the domain.
    
    if nargin < 1
        xmlFilename = 'myMesh.xml';
    end
    
    % --- Step 1: Parse the FEniCS XML file ---
    [nodes, elements, cellType] = parseFEniCSXML(xmlFilename);
    
    % --- Step 2: Automatically Plot Based on Detected Element Type ---
    figure(8); clf;
    
    switch lower(cellType)
        case 'tetrahedron'
            % Extract outer boundary triangular faces
            TR = triangulation(double(elements), double(nodes));
            bndFaces = freeBoundary(TR);
            
            trisurf(bndFaces, nodes(:,1), nodes(:,2), nodes(:,3), ...
                    'FaceColor', [0.3 0.8 0.4], 'EdgeColor', [0.2 0.2 0.2], ...
                    'FaceAlpha', 0.9, 'SpecularStrength', 0);
            title('Extruded Perimeter Mesh (Tetrahedral)');
            
        case 'hexahedron'
            % Extract all 6 quad faces per hexahedron element
            % UFC ordering: [v0, v1, v2, v3] bottom, [v4, v5, v6, v7] top
            allFaces = [
                elements(:, [1 2 3 4]);  % Bottom face (-z)
                elements(:, [5 6 7 8]);  % Top face (+z)
                elements(:, [1 2 6 5]);  % Front face (-y)
                elements(:, [2 3 7 6]);  % Right face (+x)
                elements(:, [3 4 8 7]);  % Back face (+y)
                elements(:, [4 1 5 8])   % Left face (-x)
            ];
            
            % Identify unique boundary faces (faces shared by only 1 element)
            sortedFaces = sort(allFaces, 2);
            [~, ia, ic] = unique(sortedFaces, 'rows');
            counts = accumarray(ic, 1);
            bndFaces = allFaces(ia(counts == 1), :);
            
            patch('Faces', bndFaces, 'Vertices', nodes, ...
                  'FaceColor', [0.3 0.8 0.4], 'EdgeColor', [0.2 0.2 0.2], ...
                  'FaceAlpha', 0.9, 'SpecularStrength', 0);
            title('Extruded Perimeter Mesh (Hexahedral)');
            
        case 'triangle'
            trisurf(double(elements), double(nodes(:,1)), double(nodes(:,2)), double(nodes(:,3)), ...
                    'FaceColor', [0.3 0.7 0.9], 'EdgeColor', [0.1 0.1 0.1], ...
                    'SpecularStrength', 0);
            title('Surface Triangle Mesh');
            
        otherwise
            error('Unsupported or unknown cell type: %s', cellType);
    end
    
    axis equal; grid on; view(0, 90);
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
    camlight; lighting flat;
end

function [nodes, elements, cellType] = parseFEniCSXML(xmlFilename)
    % Optimized FEniCS/DOLFIN XML parser with automatic cell-type detection
    fid = fopen(xmlFilename, 'r');
    if fid == -1
        error('Could not open file: %s', xmlFilename);
    end
    
    cellType = '';
    nodes = [];
    elements = [];
    
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line)
            break; 
        end
        line = strtrim(line);
        
        % 1. Detect Celltype
        if strncmp(line, '<mesh', 5)
            tokens = regexp(line, 'celltype="(\w+)"', 'tokens');
            if ~isempty(tokens)
                cellType = tokens{1}{1};
            end
            
        % 2. Parse Vertices
        elseif strncmp(line, '<vertices', 9)
            tokens = regexp(line, 'size="(\d+)"', 'tokens');
            numVerts = str2double(tokens{1}{1});
            
            % Matches: <vertex index="0" x="1.23" y="4.56" z="7.89" />
            formatSpec = ' <vertex index="%*d" x="%f" y="%f" z="%f" />';
            data = textscan(fid, formatSpec, numVerts);
            nodes = double([data{1}, data{2}, data{3}]);
            
        % 3. Parse Cells based on detected type
        elseif strncmp(line, '<cells', 6)
            tokens = regexp(line, 'size="(\d+)"', 'tokens');
            numCells = str2double(tokens{1}{1});
            
            switch lower(cellType)
                case 'tetrahedron'
                    formatSpec = ' <tetrahedron index="%*d" v0="%d" v1="%d" v2="%d" v3="%d" />';
                    data = textscan(fid, formatSpec, numCells);
                    elements = double([data{1}, data{2}, data{3}, data{4}] + 1);
                    
                case 'hexahedron'
                    formatSpec = ' <hexahedron index="%*d" v0="%d" v1="%d" v2="%d" v3="%d" v4="%d" v5="%d" v6="%d" v7="%d" />';
                    data = textscan(fid, formatSpec, numCells);
                    elements = double([data{1}, data{2}, data{3}, data{4}, ...
                                       data{5}, data{6}, data{7}, data{8}] + 1);
                                   
                case 'triangle'
                    formatSpec = ' <triangle index="%*d" v0="%d" v1="%d" v2="%d" />';
                    data = textscan(fid, formatSpec, numCells);
                    elements = double([data{1}, data{2}, data{3}] + 1);
                    
                otherwise
                    fclose(fid);
                    error('Unsupported cell type: %s', cellType);
            end
        end
    end
    fclose(fid);
end
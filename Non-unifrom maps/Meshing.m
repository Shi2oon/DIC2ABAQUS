function [X, Y, E11, E22, E12, Nodes, Elements, mesh] = Meshing(DATAin)
% MESHING
% Generates a quadrilateral mesh from non-uniform or incomplete DIC maps.
%
% Expected columns:
%   1: x-coordinate
%   2: y-coordinate
%   3: E11
%   4: E22
%   5: E12, optional
%
% An element is created only when:
%   1. All four corner nodes exist.
%   2. All required nodal DIC quantities are finite.
%
% Element labels are always consecutive from 1 to number of elements.

    if nargin ~= 1
        error('Meshing requires one input array.');
    end

    if size(DATAin, 2) < 4
        error(['DATAin must contain at least four columns: ', ...
               'x, y, E11 and E22.']);
    end

    if any(~isfinite(DATAin(:,1:2)), 'all')
        error('The x and y coordinates must be finite.');
    end

    % Sort by x first, then y, matching the original implementation.
    DATAin = sortrows(DATAin, [1, 2]);

    xValues = unique(DATAin(:,1), 'sorted');
    yValues = unique(DATAin(:,2), 'sorted');

    nX = numel(xValues);
    nY = numel(yValues);

    if nX < 2 || nY < 2
        error('At least two distinct x and y coordinates are required.');
    end

    % Associate each node with its position in the Cartesian coordinate map.
    [~, xIndex] = ismember(DATAin(:,1), xValues);
    [~, yIndex] = ismember(DATAin(:,2), yValues);

    nodeMap = zeros(nY, nX);

    for nodeID = 1:size(DATAin,1)

        row = yIndex(nodeID);
        col = xIndex(nodeID);

        if nodeMap(row,col) ~= 0
            error(['Duplicate DIC coordinate detected at x = %.12g, ', ...
                   'y = %.12g.'], ...
                   DATAin(nodeID,1), DATAin(nodeID,2));
        end

        nodeMap(row,col) = nodeID;
    end

    % A node is suitable for meshing only when its required DIC data exist.
    requiredData = DATAin(:,3:4);

    if size(DATAin,2) >= 5
        requiredData = DATAin(:,3:5);
    end

    validNode = all(isfinite(requiredData), 2);

    maximumNumberOfElements = (nX - 1) * (nY - 1);
    connectivity = zeros(maximumNumberOfElements, 4);

    elementCounter = 0;

    for ix = 1:(nX - 1)

        for iy = 1:(nY - 1)

            % Preserve the node ordering used by the original Meshing.m:
            %
            % node 4 -------- node 3
            %   |                |
            %   |                |
            % node 1 -------- node 2
            %
            % In coordinate terms:
            % 1: right-bottom
            % 2: right-top
            % 3: left-top
            % 4: left-bottom

            elementNodes = [ ...
                nodeMap(iy,     ix + 1), ...
                nodeMap(iy + 1, ix + 1), ...
                nodeMap(iy + 1, ix), ...
                nodeMap(iy,     ix)];

            % Do not create an element if one or more coordinate nodes
            % are absent from the DIC map.
            if any(elementNodes == 0)
                continue;
            end

            % Do not create an element if any corner contains invalid data.
            if any(~validNode(elementNodes))
                continue;
            end

            elementCounter = elementCounter + 1;
            connectivity(elementCounter,:) = elementNodes;
        end
    end

    connectivity = connectivity(1:elementCounter,:);

    if isempty(connectivity)
        error(['No valid quadrilateral elements were generated. ', ...
               'Check the DIC coordinates and missing-data mask.']);
    end

    nodeLabels = (1:size(DATAin,1)).';
    elementLabels = (1:elementCounter).';

    Nodes = [nodeLabels, DATAin(:,1:2)];
    Elements = [elementLabels, connectivity];

    % Arrange element quantities as:
    % four local nodes by number of elements.
    X = reshape(DATAin(connectivity(:),1), size(connectivity)).';
    Y = reshape(DATAin(connectivity(:),2), size(connectivity)).';

    E11 = reshape(DATAin(connectivity(:),3), size(connectivity)).';
    E22 = reshape(DATAin(connectivity(:),4), size(connectivity)).';

    if size(DATAin,2) >= 5
        E12 = reshape(DATAin(connectivity(:),5), size(connectivity)).';
    else
        E12 = [];
    end

    % Preserve the original misspelled field name for compatibility.
    mesh.winodow = [nX, nY];

    % Also provide a correctly spelled version for future use.
    mesh.window = [nX, nY];

    mesh.xy = DATAin(:,1:2).';
    mesh.Data = DATAin(:,3:end).';
    mesh.nodeMap = nodeMap;
    mesh.validNode = validNode;
    mesh.connectivity = connectivity;
end
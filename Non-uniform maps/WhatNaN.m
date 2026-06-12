function [Elements, E11, E22, E12, X, Y, answer] = ...
    WhatNaN(Elements, E11, E22, E12, X, Y)
% WHATNAN
% Removes elements containing non-finite nodal data and renumbers the
% retained elements consecutively.
%
% Supports:
%   2D quadrilateral elements:
%       Elements = [element ID, node 1, node 2, node 3, node 4]
%
%   3D brick elements:
%       Elements = [element ID, node 1, ..., node 8]
%
% Zero is treated as a valid physical value. Only NaN and Inf values cause
% an element to be removed.

    %% Basic validation

    if isempty(Elements)
        error('Elements is empty.');
    end

    numberOfElements = size(Elements, 1);
    numberOfElementNodes = size(Elements, 2) - 1;

    if ~ismember(numberOfElementNodes, [4, 8])
        error(['Unsupported element connectivity. Elements has %d columns. ', ...
               'Expected 5 columns for a 2D quadrilateral or 9 columns ', ...
               'for a 3D brick.'], size(Elements, 2));
    end

    validateElementArray( ...
        E11, numberOfElements, numberOfElementNodes, 'E11');

    validateElementArray( ...
        E22, numberOfElements, numberOfElementNodes, 'E22');

    validateElementArray( ...
        X, numberOfElements, numberOfElementNodes, 'X');

    validateElementArray( ...
        Y, numberOfElements, numberOfElementNodes, 'Y');

    if ~isempty(E12)
        validateElementArray( ...
            E12, numberOfElements, numberOfElementNodes, 'E12');
    end

    %% Find invalid elements

    invalidElement = ...
        any(~isfinite(E11), 1) | ...
        any(~isfinite(E22), 1) | ...
        any(~isfinite(X), 1)   | ...
        any(~isfinite(Y), 1);

    if ~isempty(E12)
        invalidElement = ...
            invalidElement | any(~isfinite(E12), 1);
    end

    % Ensure logical row vector.
    invalidElement = reshape(invalidElement, 1, []);

    if numel(invalidElement) ~= numberOfElements
        error(['The generated invalid-element mask has the wrong size. ', ...
               'Expected %d values but obtained %d.'], ...
               numberOfElements, numel(invalidElement));
    end

    keepElement = ~invalidElement;

    %% Remove invalid elements and corresponding element data

    Elements = Elements(keepElement, :);

    E11 = E11(:, keepElement);
    E22 = E22(:, keepElement);

    X = X(:, keepElement);
    Y = Y(:, keepElement);

    if ~isempty(E12)
        E12 = E12(:, keepElement);
    end

    if isempty(Elements)
        error(['All elements were removed because every element contained ', ...
               'at least one NaN or Inf value.']);
    end

    %% Renumber the retained elements

    % The first column contains the Abaqus element label.
    Elements(:, 1) = (1:size(Elements, 1)).';

    %% Final consistency checks

    retainedElementCount = size(Elements, 1);

    assert(size(E11, 2) == retainedElementCount, ...
        'E11 is inconsistent with the retained element table.');

    assert(size(E22, 2) == retainedElementCount, ...
        'E22 is inconsistent with the retained element table.');

    assert(size(X, 2) == retainedElementCount, ...
        'X is inconsistent with the retained element table.');

    assert(size(Y, 2) == retainedElementCount, ...
        'Y is inconsistent with the retained element table.');

    if ~isempty(E12)
        assert(size(E12, 2) == retainedElementCount, ...
            'E12 is inconsistent with the retained element table.');
    end

    expectedElementLabels = (1:retainedElementCount).';

    assert(isequal(Elements(:, 1), expectedElementLabels), ...
        'The retained elements were not renumbered correctly.');

    %% Report result

    numberRemoved = sum(invalidElement);

    if numberRemoved > 0
        answer = 'Y';

        fprintf(['WhatNaN removed %d invalid elements and retained %d ', ...
                 '%d-node elements.\n'], ...
                 numberRemoved, retainedElementCount, ...
                 numberOfElementNodes);
    else
        answer = 'N';

        fprintf(['WhatNaN found no invalid elements. Retained %d ', ...
                 '%d-node elements.\n'], ...
                 retainedElementCount, numberOfElementNodes);
    end
end


function validateElementArray( ...
    array, numberOfElements, numberOfElementNodes, arrayName)
% VALIDATEELEMENTARRAY
% Checks that an element-associated array has one column per element and
% one row per local element node.

    if isempty(array)
        error('%s is empty.', arrayName);
    end

    if size(array, 2) ~= numberOfElements
        error(['%s has %d columns, but Elements contains %d rows. ', ...
               'There must be one column per element.'], ...
               arrayName, size(array, 2), numberOfElements);
    end

    if size(array, 1) ~= numberOfElementNodes
        error(['%s has %d rows, but the element connectivity contains ', ...
               '%d nodes per element.'], ...
               arrayName, size(array, 1), numberOfElementNodes);
    end
end

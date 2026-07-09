function [Elements, E11, E22, E12, X, Y, answer] = ...
    WhatNaN(Elements, E11, E22, E12, X, Y)
% WHATNAN
% Removes finite elements that contain invalid E11 or E22 values at any
% local element node.
%
% The function supports:
%   CPS4/CPE4 elements:
%       Elements = [elementID, n1, n2, n3, n4]
%
%   C3D8 elements:
%       Elements = [elementID, n1, n2, n3, n4, n5, n6, n7, n8]
%
% The strain and coordinate arrays must have the form:
%
%       numberOfLocalNodes x numberOfElements
%
% Geometry removal is based only on E11 and E22. Zero strain is considered
% valid and statistical outlier detection is not performed.

    %% Validate the element table

    if isempty(Elements)
        error('WhatNaN:EmptyElements', ...
              'The element connectivity table is empty.');
    end

    numberOfElements = size(Elements, 1);
    numberOfLocalNodes = size(Elements, 2) - 1;

    if ~ismember(numberOfLocalNodes, [4, 8])
        error('WhatNaN:UnsupportedElementType', ...
             ['Unsupported element connectivity with %d local nodes. ', ...
              'Expected four-node quadrilaterals or eight-node bricks.'], ...
              numberOfLocalNodes);
    end

    %% Validate element-associated arrays

    validateElementArray( ...
        E11, numberOfLocalNodes, numberOfElements, 'E11');

    validateElementArray( ...
        E22, numberOfLocalNodes, numberOfElements, 'E22');

    validateElementArray( ...
        X, numberOfLocalNodes, numberOfElements, 'X');

    validateElementArray( ...
        Y, numberOfLocalNodes, numberOfElements, 'Y');

    if ~isempty(E12)
        validateElementArray( ...
            E12, numberOfLocalNodes, numberOfElements, 'E12');
    end

    %% Identify invalid strain points

    invalidE11Point = isnan(E11) | isinf(E11);
    invalidE22Point = isnan(E22) | isinf(E22);

    % A point is invalid when either normal strain component is invalid.
    invalidStrainPoint = invalidE11Point | invalidE22Point;

    %% Convert the nodal mask to an element mask

    % Each column corresponds to one element.
    % Remove the element when any local node is invalid.
    invalidElement = any(invalidStrainPoint, 1);

    % Force the mask to be a logical row vector.
    invalidElement = logical(reshape(invalidElement, 1, []));

    if numel(invalidElement) ~= numberOfElements
        error('WhatNaN:InvalidMaskSize', ...
             ['The invalid-element mask contains %d entries, but the ', ...
              'element table contains %d elements.'], ...
              numel(invalidElement), numberOfElements);
    end

    keepElement = ~invalidElement;

    oldElementLabels = Elements(:, 1);
    removedElementLabels = oldElementLabels(invalidElement);

    %% Remove invalid elements and corresponding data

    Elements = Elements(keepElement, :);

    E11 = E11(:, keepElement);
    E22 = E22(:, keepElement);

    X = X(:, keepElement);
    Y = Y(:, keepElement);

    if ~isempty(E12)
        E12 = E12(:, keepElement);
    end

    if isempty(Elements)
        error('WhatNaN:AllElementsRemoved', ...
             ['Every element contained at least one invalid E11 or E22 ', ...
              'point. Check the strain data and coordinate ordering.']);
    end

    %% Renumber retained elements consecutively

    % Preserve node labels. Only element labels are changed.
    Elements(:, 1) = (1:size(Elements, 1)).';

    %% Final consistency checks

    retainedElementCount = size(Elements, 1);

    assert(size(E11, 2) == retainedElementCount, ...
        'E11 no longer corresponds to the element table.');

    assert(size(E22, 2) == retainedElementCount, ...
        'E22 no longer corresponds to the element table.');

    assert(size(X, 2) == retainedElementCount, ...
        'X no longer corresponds to the element table.');

    assert(size(Y, 2) == retainedElementCount, ...
        'Y no longer corresponds to the element table.');

    if ~isempty(E12)
        assert(size(E12, 2) == retainedElementCount, ...
            'E12 no longer corresponds to the element table.');
    end

    expectedLabels = (1:retainedElementCount).';

    assert(isequal(Elements(:, 1), expectedLabels), ...
        'Element renumbering failed.');

    assert(~any(isnan(E11(:)) | isinf(E11(:))), ...
        'Invalid values remain in E11 after element removal.');

    assert(~any(isnan(E22(:)) | isinf(E22(:))), ...
        'Invalid values remain in E22 after element removal.');

    %% Report removal

    numberRemoved = sum(invalidElement);

    if numberRemoved > 0
        answer = 'Y';

        fprintf('\n');
        fprintf('WhatNaN geometry trimming:\n');
        fprintf('  Original elements: %d\n', numberOfElements);
        fprintf('  Removed elements:  %d\n', numberRemoved);
        fprintf('  Retained elements: %d\n', retainedElementCount);
        fprintf('  Local nodes/element: %d\n', numberOfLocalNodes);

        fprintf('  Original labels of removed elements:\n');
        fprintf('  ');

        fprintf('%d ', removedElementLabels);

        fprintf('\n\n');
    else
        answer = 'N';

        fprintf('\n');
        fprintf(['WhatNaN found no NaN or Inf values in E11 or E22. ', ...
                 'All %d elements were retained.\n\n'], ...
                 numberOfElements);
    end
end


function validateElementArray( ...
    array, expectedRows, expectedColumns, arrayName)
% VALIDATEELEMENTARRAY
% Confirms that an element-associated array has one column per element and
% one row per local element node.

    if isempty(array)
        error('WhatNaN:EmptyArray', ...
              '%s is empty.', arrayName);
    end

    if size(array, 1) ~= expectedRows
        error('WhatNaN:InvalidRowCount', ...
             ['%s contains %d rows. It must contain %d rows, one for ', ...
              'each local element node.'], ...
              arrayName, size(array, 1), expectedRows);
    end

    if size(array, 2) ~= expectedColumns
        error('WhatNaN:InvalidColumnCount', ...
             ['%s contains %d columns, while Elements contains %d rows. ', ...
              'There must be one column per element.'], ...
              arrayName, size(array, 2), expectedColumns);
    end
end
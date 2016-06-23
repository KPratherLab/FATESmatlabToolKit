function [lia,locb] = ismemberR2012a_Get_spectrum_peak_script(a,b,options)
%this is a copy of ismemberR2012a (a built in matlab funciton)
%that has been altered to remove the 
%uniqueness checks for inputs a and b.  This is a way to speed up
%get_spectrum_peak_script as the PID inputs (a and b) have already been
%checked for uniqueness externally.  The uniqueness check is very time
%consuming for long lists and ends up increasing the get_spectrum call 
%non linearly with length of a. The only change involved commenting out lines 94-96 and 98-103 
%and adding two new lines below it creating the iaC and ib variables that
%are called later on.

%ismemberBuiltInTypes also had to be copied in from the built in ismember
%MATLAB function

% 'R2012a' flag implementation

% Error check flag
if nargin == 2
    byrow = false;
else
    byrow = options > 0;
end

doBuiltinTypes = true;
% Check that one of A and B is double if A and B are non-homogeneous. Do a
% separate check if A is a heterogeneous object and only allow a B
% that is of the same root class.
if ~(isa(a,'handle.handle') || isa(b,'handle.handle'))
    if ~strcmpi(class(a),class(b))
        if isa(a,'matlab.mixin.Heterogeneous') && isa(b,'matlab.mixin.Heterogeneous')
            rootClassA = meta.internal.findHeterogeneousRootClass(a);
            if isempty(rootClassA) || ~isa(b,rootClassA.Name)
                error(message('MATLAB:ISMEMBER:InvalidInputsDataType',class(a),class(b)));
            end
            doBuiltinTypes = false;
        elseif ~(strcmpi(class(a),'double') || strcmpi(class(b),'double'))
            error(message('MATLAB:ISMEMBER:InvalidInputsDataType',class(a),class(b)));
        end
    end
end

if ~byrow
    if ~(isa(a,'opaque') || isa(b,'opaque')) && doBuiltinTypes
        % Builtin types
        if nargout > 1
            [lia,locb] = ismemberBuiltinTypes(a,b);
        else
            lia = ismemberBuiltinTypes(a,b);
        end
    else
        % Handle objects
        if nargout > 1
            [lia,locb] = ismemberClassTypes(a,b);
        else
            lia = ismemberClassTypes(a,b);
        end
    end
else    % 'rows' case
    if ~(ismatrix(a) && ismatrix(b))
        error(message('MATLAB:ISMEMBER:NotAMatrix'));
    end
    
    [rowsA,colsA] = size(a);
    [rowsB,colsB] = size(b);
    
    % Automatically pad strings with spaces
    if ischar(a) && ischar(b),
        b = [b repmat(' ',rowsB,colsA-colsB)];
        a = [a repmat(' ',rowsA,colsB-colsA)];
    elseif colsA ~= colsB
        error(message('MATLAB:ISMEMBER:AandBColnumAgree'));
    end
    
    % Empty check for 'rows'.
    if rowsA == 0 || rowsB == 0
        lia = false(rowsA,1);
        locb = zeros(rowsA,1);
        return
    end
    
    % General handling for 'rows'.
    
    % Duplicates within the sets are eliminated
    if (rowsA == 1)
        uA = repmat(a,rowsB,1);
        d = uA(1:end,:)==b(1:end,:);
        d = all(d,2);
        lia = any(d);
        if nargout > 1
            if lia
                locb = find(d, 1, 'first');
            else
                locb = 0;
            end
        end
        return;
%     else %These lines are eliminated since tmpPID does not need to be unique
%     . This speeds up get_spectrum_peak_script 
%         [uA,~,icA] = unique(a,'rows','sorted');
    end
%     if nargout <= 1 %These lines are eliminated since unique is performed outside
%     this function already for UKID. This speeds up get_spectrum_peak_script
%         uB = unique(b,'rows','sorted');
%     else
%         [uB,ib] = unique(b,'rows','sorted');
%     end
    icA = 1:size(a,1); %new line 
    ib = 1:size(b,1); %new line
    % Sort the unique elements of A and B, duplicate entries are adjacent
    [sortuAuB,IndSortuAuB] = sortrows([a;b]);
    
    % Find matching entries
    d = sortuAuB(1:end-1,:)==sortuAuB(2:end,:);     % d indicates matching entries
    d = all(d,2);                 
                   % Finds the index of matching entries
    ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C
    
    if nargout <= 1
        lia = ismemberBuiltinTypes(icA,ndx1);           % Find repeats among original list
    else
        szuA = size(a,1);
        [lia,locb] = ismemberBuiltinTypes(icA,ndx1);    % Find locb by using given indices
        d = find(d);
        newd = d(locb(lia));                    % NEWD is D for non-unique A
        where = ib(IndSortuAuB(newd+1)-szuA);   % Index values of uB through UNIQUE
        locb(lia) = where;                      % Return first or last occurrence of A within B
    end
end
end

function [lia,locb] = ismemberBuiltinTypes(a,b)
% General handling.
% Use FIND method for very small sizes of the input vector to avoid SORT.
if nargout > 1
    locb = zeros(size(a));
end
% Handle empty arrays and scalars.  
numelA = numel(a);
numelB = numel(b);
if numelA == 0 || numelB <= 1
    if numelA > 0 && numelB == 1
        lia = (a == b);
        if nargout > 1
            % Use DOUBLE to convert logical "1" index to double "1" index.
            locb = double(lia);
        end
    else
        lia = false(size(a));
    end
    return
end

scalarcut = 5;
if numelA <= scalarcut
    lia = false(size(a));
    if nargout <= 1
        for i=1:numelA
            lia(i) = any(a(i)==b(:));   % ANY returns logical.
        end
    else
        for i=1:numelA
            found = a(i)==b(:);  % FIND returns indices for LOCB.
            if any(found)
                lia(i) = true;
                found = find(found);
                locb(i) = found(1);
            end
        end
    end
else
    % Use method which sorts list, then performs binary search.
    % Convert to full to work in C helper.
    if issparse(a)
        a = full(a);
    end
    if issparse(b)
        b = full(b);
    end
    
    if (isreal(b))
        % Find out whether list is presorted before sort
        % If the list is short enough, SORT will be faster than ISSORTED
        % If the list is longer, ISSORTED can potentially save time
        checksortcut = 1000;
        if numelB > checksortcut
            sortedlist = issorted(b(:));
        else
            sortedlist = 0;
        end
        if nargout > 1
            if ~sortedlist
                [b,idx] = sort(b(:));
            end
        elseif ~sortedlist
            b = sort(b(:));
        end
    else
        sortedlist = 0;
        [~,idx] = sort(real(b(:)));
        b = b(idx);
    end
    
    % Use builtin helper function ISMEMBERHELPER:
    % [LIA,LOCB] = ISMEMBERHELPER(A,B) Returns logical array LIA indicating
    % which elements of A occur in B and a double array LOCB with the
    % locations of the elements of A occuring in B. If multiple instances
    % occur, the first occurence is returned. B must be already sorted.
    
    if ~isobject(a) && ~isobject(b) && isnumeric(a)
        if (isnan(b(end)))
            % If NaNs detected, remove NaNs from B.
            b = b(~isnan(b(:)));
        end
        if nargout <= 1
            lia = builtin('_ismemberhelper',a,b);
        else
            [lia, locb] = builtin('_ismemberhelper',a,b);
        end
    else %(a,b, are some other class like gpuArray, syb object)
        lia = false(size(a));
        if nargout <= 1
            for i=1:numelA
                lia(i) = any(a(i)==b(:));   % ANY returns logical.
            end
        else
            for i=1:numelA
                found = a(i)==b(:); % FIND returns indices for LOCB.
                if any(found)
                    lia(i) = true;
                    found = find(found);
                    locb(i) = found(1);
                end
            end
        end
    end
    if nargout > 1 && ~sortedlist
        % Re-reference locb to original list if it was unsorted
        locb(lia) = idx(locb(lia));
    end
end
end
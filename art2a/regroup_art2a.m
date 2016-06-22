function [outPartID, outWeightMatrix, RegroupTable ] = regroup_art2a(inPID, inWeightMatrix, RegroupFactor)
%Regroup uses the dot product value and compares it the threshold value set by the RegroupFactor to 
%combine matching clusters.  After regrouping, the clusters are sorted.
%
%Call as [outPartID, outWeightMatrix, RegroupTable] = Regroup(inPID, inWeightMatrix, RegroupFactor); where
%		inWeightMatrix is the characteristic weight matrix for the particles in inPID to be compared and regrouped.
%       Column n of inWeightMatrix corresponds to the particles in cell n of 
%       inPID. 
%		inPID is the ART2a result corresponding to inWeightMatrix.  
%       inPID is a cell array where each cell contains a nx2 matrix of
%       particle identifiers.
%       inPID{x}(:,1) = InstID
%       inPID{x}(:,2) = PartID
%       RegroupFactor is the threshold used for regrouping.
%       
%       outPartID is the regrouped ART2a result.
%       outPartID is a cell array where each cell contains of nx2 matrix of
%       particle identifiers
%       outWeightMatrix is the regrouped weight matrix.
%       RegroupTable gives a list of the original clusters that were combined and the factor at which they matched.
%
% Shields, Sodeman, Wenzel 28 Oct 2003.
% Slight Modification made on 28 Jan 2004.
% Minor fix of FirstPart on Feb 18, 2004.

% Changed input and output Weight/Area Matrixes so that each cluster WM
% is stored in a column rather than a row.  This decreases the run time.
% Also made various coding changes to improve run time. 
% Camille Sultana May 2014
% Modified to be compatible with FATES Toolkit
% Camille Sultana 2016

fprintf('INFO, regroup, \n');

%% check input variables
if nargin < 3
   error('Wrong number of input arguments');
end
if ~isreal(inWeightMatrix);
   error('Inputs must be matrices of real numbers');
end

if isempty(inWeightMatrix);
   error('Input matrix must contain data.');
end
if ~exist('RegroupFactor', 'var');
    RegroupFactor = 0
else
    if RegroupFactor < 0 | RegroupFactor > 1.0;
        error('RegroupFactor must be a number between 0 and 1.0');
    end
end  

%% 
% set up variables
NoPart = []; 
K = 0;
Table2{1} = 'I';
Table2{2} = 'J';
Table2{3} = 'CompareValue';
Table2 = Table2';

%do first comparison of weight matrices
Comparison = inWeightMatrix' * inWeightMatrix; %do first comparison
Comparison(logical(eye(size(Comparison)))) = 0; %set diagonal to 0
[BigValue, indxBIG] = max(Comparison(:)); %find the max value in comparison
[I, J] = ind2sub(size(Comparison),indxBIG); %find the location of the max value in comparison

%combine clusters if dot product is > threshhold
while BigValue > RegroupFactor;
    K = K + 1; % counter #. how many times a cluster gets regrouped.
    clear outPID
    %edit weight matrix
    inWeightMatrix1b = (inWeightMatrix(:,I)*size(inPID{I},1) + inWeightMatrix(:,J)*size(inPID{J},1))/(size(inPID{I},1) + size(inPID{J},1)); %combine weight matrices
    inWeightMatrix1b = inWeightMatrix1b./norm(inWeightMatrix1b); %normalize
    inWeightMatrix(:,I) = inWeightMatrix1b(:,1); %add new weight vector to weight matrix
    inWeightMatrix(:,J) = 0; %set second cluster weight vector to 0
    
    %edit PID
    [outPID] = [inPID{I}; inPID{J}]; %combine particles in matched clusters
    inPID{I} = outPID;
    inPID{J} = NoPart; %clear pids second cluster that got combined
    
    %keep record of clusters that were combined
    Table{1} = I;
    Table{2} = J;
    Table{3} = BigValue;
    Table1 = Table';
    Table2 = [Table2, Table1];
    
    %do comparison of weight matrices again
    Comparison          = inWeightMatrix' * inWeightMatrix;
    Comparison(logical(eye(size(Comparison)))) = 0;
    [BigValue, indxBIG] = max(Comparison(:));
    [I, J]              = ind2sub(size(Comparison),indxBIG);
end

RegroupTable = Table2';
fprintf('INFO,regroup, %i merging happened\n',K);

%% sort results 
noPID             = cell2mat(cellfun(@isempty, inPID, 'uniformoutput', false)); %remove empty clusters
inPID = inPID(~noPID);
inWeightMatrix = inWeightMatrix(:,~noPID);
TempLength = cell2mat(cellfun(@(x) size(x,1), inPID, 'uniformoutput', false));
[partIDCount, SORTid] = sort(TempLength,'descend'); %sort by length
outPartID = inPID(SORTid); %reorder clusters getting rid of empty ones
for i = 1:length(outPartID)
    outPartID{i} = sortrows(outPartID{i});
end
outWeightMatrix = inWeightMatrix(:,SORTid); %reorder weight matrix

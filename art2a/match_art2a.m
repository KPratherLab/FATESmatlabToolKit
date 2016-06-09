function [outPartID, Final_Part2Nrn_Idx, Final_Part2Nrn_Prox] = match_art2a(inPartID, WeightMatrix, Polarity, VigilanceFactor, Exclusive)
% Call as [outPartID] = match_art2a_CS(inPartID, WeightMatrix, Polarity, VigilanceFactor, Exclusive)
% where outPartID is a cell array where each cell contains the PartIDs the particles that are within the 
%          specified VigilanceFactor to the corresponding cluster in the WeightMatrix, such that column M of 
%          WeightMatrix corresponds to the	particles in cell M of outPartID.
%
%       inPartID is a nx2 matrix of particle identifiers.
%       inPartID(:,1) = InstID
%       inPartID(:,2) = PartID
%       WeightMatrix is the set of cluster weight vectors (mass spectra) to match particles to.
%       Polarity is a switch to decide if there polarity should be positive, negative or dual ion.  
%           If Polarity is set to 0, negative spectra are mandated, 1, positive spectra are mandated, 
%           and 2, both polarities are considered.  Although dual polarity is the default.
%       VigilanceFactor is a double, recommended to be 0.7 but may be between 0 and 1, noninclusive
%		Exclusive (optional), if set to 1, puts each particle PartID into the most closely matched cluster 
%           in the WeigthMatrix and does not match the same particle with multiple clusters in the WeightMatrix.
%           If the particle matches equally well with 2 or more clusters in the WeightMatrix, the first one is
%           returned - so that the order of the clusters in the WeightMatrix will have some effects on the 
%           results of this function when the Exclusive option is selected.
%		  

%Copyright 1999-2001, David Fergenson and Sylvia Pastor.  All rights reserved by the Regents of the
%University of California.  Although this software is compatible with the 
%YAADA data analysis toolkit, which is governed by the Gnu Public License,
%this program, itself, is not freely distributable.
% S.H.Pastor  21 March 01 Modified for Yaada0.93d   
% ***EQUAL WEIGHT GIVEN TO EACH POLARITY, 25April2001 modified by DF***

%  Check for errors in inputs

% Modified by Alberto Cazorla. It runs the division in parts in this
% function and call to the original (subfunction below this one) renamed as
% martch_art2a_part
% November 8th, 2010

fprintf('INFO, match art2a, starting\n');

% Check for errors in inputs
if nargin < 2
    error('Too few input arguments.');
elseif nargin > 5
    error('Too many input arguments.');
end

if nargout > 3
    error('Too many output arguments.');
end

if ~(size(inPartID,2)==2)
    whos inPartID
    error('inPartID must be a Nx2 matrix of InstIDs,PartIDs');
end

if exist('WeightMatrix','var');
    if ~isreal(WeightMatrix);
        error('inWeightMatrix is non-numerical');
    end
end

%  Set variables if not specified
if ~exist('VigilanceFactor', 'var')
    VigilanceFactor = 0.7;
end

if ~exist('Exclusive', 'var')
   Exclusive = 1;
end

if ~exist('Polarity', 'var');
    Polarity = 2;
else
    if Polarity > 2 | Polarity < 0 | ~isreal(Polarity);
        error('Polarity must be a number between 0 and 2.  See help file');
    end
end

%%
%this is set up to only do 100k at a time so that getintspectrum
%only loads 100k part-peak spectrums at a time
%set variables
Final_Part2Nrn_Idx =[];
Final_Part2Nrn_Prox=[];
Part_setsize=100000;  %try 100k

%loop through matching 100K particles at a time
if length (inPartID) < Part_setsize %only do once if < 50K particles
    [outPartID Final_Part2Nrn_Idx Final_Part2Nrn_Prox]= match_art2a_part(inPartID, WeightMatrix, Polarity, VigilanceFactor, Exclusive);
    disp(5)
    whos Final*
else
    numInt = ceil((length (inPartID)) / Part_setsize); %find number of 100K chunks
    for i = 1:numInt %get 100K lists of particles
        if i < numInt
            inPID2 = inPartID((i-1)*Part_setsize + 1:i*Part_setsize,:);            
        else
            inPID2 = inPartID((1+Part_setsize*(i-1)):length(inPartID),:);            
        end
        [matched Part2Nrn_Idx Part2Nrn_Prox]= match_art2a_part(inPID2, WeightMatrix, Polarity, VigilanceFactor, Exclusive); %perform matching
           %whos *Part2Nrn*
        if i == 1
           outPartID          = matched;
           if (size(Part2Nrn_Idx,1)==1)
                Final_Part2Nrn_Idx =Part2Nrn_Idx';
           else  Final_Part2Nrn_Idx =Part2Nrn_Idx;
           end;
           Final_Part2Nrn_Prox=Part2Nrn_Prox;
        else
           outPartID          = cellfun(@(X,Y) [X; Y], matched, outPartID, 'UniformOutput',false);
           if (size(Part2Nrn_Idx,1)==1)
                    Final_Part2Nrn_Idx =[Final_Part2Nrn_Idx;  Part2Nrn_Idx'];
           else     Final_Part2Nrn_Idx =[Final_Part2Nrn_Idx;  Part2Nrn_Idx];
           end;
           Final_Part2Nrn_Prox=[Final_Part2Nrn_Prox; Part2Nrn_Prox];
        end
    end
end

for i = 1:length(outPartID)
    outPartID{i} = sortrows(outPartID{i});
end

function [outPartID, Part2Nrn_Idx, Part2Nrn_Prox] = match_art2a_part(inPartID, WeightMatrix, Polarity, VigilanceFactor, Exclusive)
% % Call as [outPartID] = match_art2a_part_CS(inPartID, WeightMatrix, Polarity, VigilanceFactor, Exclusive)
% where outPartID is a cell array where each cell contains the PartIDs the particles that are within the
%          specified VigilanceFactor to the corresponding cluster in the WeightMatrix, such that row M of
%          WeightMatrix corresponds to the	particles in cell M of outPartID.
%		  inPartID is an object of PartIDs to be considered in the analysis
%       WeightMatrix is the set of clusters to match particles to.
%       Polarity is a switch to decide if there polarity should be positive, negative or dual ion.
%           If Polarity is set to 0, negative spectra are mandated, 1, positive spectra are mandated,
%           and 2, both polarities are considered.  Although dual polarity is the default.
%       VigilanceFactor is a double, recommended to be 0.7 but may be between 0 and 1, noninclusive
%		  Exclusive (optional), if set to 1, puts each particle PartID into the most closely matched cluster
%           in the WeigthMatrix and does not match the same particle with multiple clusters in the WeightMatrix.
%           If the particle matches equally well with 2 or more clusters in the WeightMatrix, the first one is
%           returned - so that the order of the clusters in the WeightMatrix will have some effects on the
%           results of this function when the Exclusive option is selected.
%

%Copyright 1999-2001, David Fergenson and Sylvia Pastor.  All rights reserved by the Regents of the
%University of California.  Although this software is compatible with the
%YAADA data analysis toolkit, which is governed by the Gnu Public License,
%this program, itself, is not freely distributable.
% S.H.Pastor  21 March 01 Modified for Yaada0.93d
% ***EQUAL WEIGHT GIVEN TO EACH POLARITY, 25April2001 modified by DF***

%PFR added output the Par2Nrn index to say to which neuron (cluster) each particle was
%clustered


fprintf('INFO, match art2a part, starting\n');

%%
% Check for errors in inputs
if nargin < 2
    error('Too few input arguments.');
elseif nargin > 5
    error('Too many input arguments.');
end

if nargout > 3
    error('Too many output arguments.');
end

if ~(size(inPartID,2)==2)
    whos inPartID
    error('inPartID must be a Nx2 matrix of InstIDs,PartIDs');
end

%  Set Vigilance if not specified
if ~exist('VigilanceFactor', 'var')
    VigilanceFactor = 0.7;
end

if exist('WeightMatrix','var');
    if ~isreal(WeightMatrix);
        error('inWeightMatrix is non-numerical');
    end
end

%%
%set up variables
if ~exist('Polarity', 'var');
    Polarity = 2;
else
    if Polarity > 2 | Polarity < 0 | ~isreal(Polarity);
        error('Polarity must be a number between 0 and 2.  See help file');
    end
end

if ~exist('Exclusive', 'var')
    Exclusive = 1;
end

switch Polarity
    case 0;
        MaxMz = size(WeightMatrix,1);
    case 1;
        MaxMz = size(WeightMatrix,1);
    case 2;
        MaxMz = size(WeightMatrix,1)/2;
end;

NumSeeds = size(WeightMatrix, 2);

%%
%get spectra
switch Polarity
    case 0;
        [AreaMatrix] = get_int_spectrum_SUM(inPartID,MaxMz,'Area',0);
    case 1;
        [~,AreaMatrix] = get_int_spectrum_SUM(inPartID,MaxMz,'Area',1);
    case 2;
        [NegArea,PosArea] = get_int_spectrum_SUM(inPartID,MaxMz,'Area',2);
        hasdataPos = any(PosArea,1);
        hasdataNeg = any(NegArea,1);
        PosArea = normc(PosArea);
        NegArea = normc(NegArea);
        PosArea(:,~hasdataPos) = 0;
        NegArea(:,~hasdataNeg) = 0;
        AreaMatrix = [PosArea; NegArea];
end;

if isempty(AreaMatrix)
    error('PIDs did not generate AreaMatrix with any peaks');
end

%normalize spectra
hasdata = any(AreaMatrix,1);
AreaMatrix = AreaMatrix(:,hasdata);
inPartID = inPartID(hasdata,:);
AreaMatrix = normc(AreaMatrix);
NumSpectra = size(inPartID,1);

%%
%start matching

%make sure WeightMatrix and AreaMatrix have same number of m/z
Temp = size(WeightMatrix, 1) - size(AreaMatrix, 1);
if Temp > 0;
    error(['The weight matrix was too wide by ', num2str(Temp), '.']);
elseif Temp < 0;
    error(['The weight matrix was too narrow by ', num2str(Temp), '.']);
end
%Compare the AreaMatrix from inPartID to the WeightMatrix
ProximityMatrix              = AreaMatrix' * WeightMatrix; %do dot product
%find cluster that best matches each particle
[Part2Nrn_Prox, Part2Nrn_Idx]  = max(ProximityMatrix,[],2); %MaxProx is value, NumWtVect is cluster
Part2Nrn_Idx(Part2Nrn_Prox < VigilanceFactor) = 0; %if dot product is less VigFac than no match to cluster

%sort particles into clusters
if Exclusive;   %sort into best matched cluster
    for I = 1:NumSeeds;  %NumSeesds is # clusters (ie neurons)
        outPartID{I} = inPartID(Part2Nrn_Idx==I,:);
    end
else %sort into any cluster with dot product exceeding threshold
    for I = 1:NumSeeds;
        Threshold        = ProximityMatrix >= VigilanceFactor; %find every cluster where dot product exceeds threshold
        outPartID{I}     = inPartID(Threshold(:,I),:);
    end
end

disp('end of part subfnct');
disp(1)
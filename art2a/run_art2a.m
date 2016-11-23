function [outPartID, outWeightMatrix, partIDCount, trackNeuron, avgPROXALL, PDistALL] = run_art2a(inPartID, Polarity, MaxMass, LearningRate, VigilanceFactor, MaxIteration, NegArea, PosArea)
% RUN_ART2A: ART2a implementation for the YAADA Database.
% Call as function [outPartID, outWeightMatrix, partIDCount, trackNeuron, PDistALL] = 
% run_art2a_CS_pr_v2(inPartID, Polarity, MaxMass, LearningRate, VigilanceFactor, MaxIteration)
%
%   inPartID is set of paritcle identifiers stored as a nx2 matrix
%       inPartID(:,1) = InstID
%       inPartID(:,2) = PartID%
%   Polarity is a switch to decide if the polarity should be positive, negative
%       or dual ion.  If Polarity is set to 0, negative spectra are mandated,
%       1, positive spectra are mandated, and 2, both polarities are considered.
%       Dual Polarity (2) is default.
%   MaxMass (optional) is the maximum mass to charge ratio (default is
%               the YAADA.MaxMZ),
%   LearningRate (optional) is the learning rate (default is 0.05),
%   VigilanceFactor (optional)is the vigilance factor (default is 0.8),
%   MaxIteration (optional) is the maximum number of iterations (default is 20)
% 
% outPartID is a cell array where each cell contains a nx2 matrix of
%   particle identifiers.
%   outPartID{x}(:,1) = InstID
%   outPartID{x}(:,2) = PartID%
% outWeightMatrix is the characteristic weight matrix of the corresponding particles
%   such that column M of outWeightMatrix corresponds to the particles in
%   cell M of outPartID. If Dual Polarity (2) was used
%   outWeightMatrix(1:MaxMass,:) corresponds to m/z 1 to MaxMass and
%   outWeightMatrix(MaxMass:end,:) corresponsed to m/z -1 to -MaxMass.
% partIDCount is a vector with the number of particles per cluster
% avgPROXALL is a cell array that gives the average proximity of each
%   particle to the current WM of its respective cluster after each
%   iteration.    
% PDistALL is a cell array that gives the average pair wise
%   (euclidean)distance for all particles within a cluster after each
%   iteration.  SEE FOLLOWING NOTE. PDistALL increases the runtime x3
%   so only use when you are diagnosing the 'tightness' of the
%   clusters. Search 'PDistALL' to find lines to uncomment.
% trackNeuron gives the fraction of particles moved to different
%   clusters between iterations.
%
%%
%Version 1.2, Davio, January 24, 2000

%Copyright 1999-2001, David Fergenson.  All rights reserved by the Regents of the
%University of California.  Although this software is compatible with the 
%YAADA data analysis toolkit, which is governed by the Gnu Public License,
%this program, itself, is not freely distributable.

% Modified to include as parameter the learning rate, vigilance factor and
% maximum mass.
% Modified to provide the number of particles per cluster
% Also sort and truncate the outputs (no need to call truncate_art2a_result)
% Alberto Cazorla. 2010-10-8

% Rewrote script to improve run-time speed and also include metrics to track cluster
% tightness and homogeneity (avgPROXALL, PDistALL, trackNeuron).
% Removed the fixed seeds capability (use match_art2a instead).
% Changed outWeightMatrix so that each cluster WM
% is stored in a column rather than a row.
% Camille Sultana May 2014

fprintf('INFO, run_art2a_CS, starting\n');

%If you want the run to be random (to characterize ART-2a for example), then comment out the next line (rng).
rng('default');

%%
%check input variables
if nargin == 0
   error('Too few input arguments.');
elseif nargin > 8  
    error('Too many input arguments.');
end
if nargout > 6
   error('Too many output arguments.');
end

if nargout > 3 %if more than three output variables then will want to calc trackNeuron
    trackN = 1;
else
    trackN = 0;
end

if nargout > 4 %if more than 4 output variables then will want to calc avgProxall
    avgP = 1;
else
    avgP = 0;
end

if nargout > 5 %if more than 5 output variables then will want to calc PDistall
    pD = 1;
else
    pD = 0;
end

if (~(size(inPartID,2)==2))  %check for Nx2 matrix
  error('Invalid PID');
end

if ~exist('MaxMass', 'var')
   MaxMass = YAADA.MaxMZ;
end

if ~exist('LearningRate', 'var')
   LearningRate = 0.05;
end

if ~exist('VigilanceFactor', 'var')
   VigilanceFactor = 0.8;
end

if ~exist('MaxIteration', 'var')
   MaxIteration = 20;
end

if ~exist('Polarity', 'var');
   Polarity = 2;
elseif Polarity > 2 | Polarity < 0 | ~isreal(Polarity);
      error('Polarity must be a number between 0 and 2.  See help file');
end

%%
%get area table for input PID
switch Polarity
    case 0;
        if ~exist('NegArea', 'var')
            [AreaMatrix] = get_int_spectrum_SUM(inPartID,MaxMass,'Area',0); 
            fprintf('INFO, run art2a, getting int spectrum for Neg Area');
        end
       
    case 1;
        if ~exist('PosArea', 'var')
            [~,AreaMatrix] = get_int_spectrum_SUM(inPartID,MaxMass,'Area',1);
            fprintf('INFO, run art2a, getting int spectrum for Pos Area');
        end
        
    case 2;
        if (~exist('NegArea', 'var') || ~exist('PosArea','var'))
            [NegArea, PosArea] = get_int_spectrum_SUM(inPartID,MaxMass,'Area',2);
            fprintf('INFO, run art2a, getting int spectrum for Neg and Pos Area');
        end;
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

%Make sure it doesn't choke on the intersects.
hasdata = any(AreaMatrix,1);
AreaMatrix = AreaMatrix(:,hasdata);
inPartID = inPartID(hasdata,:);
AreaMatrix = normc(AreaMatrix);
NumPID = size(inPartID,1);

%%
%run art2a

% set up variables
    Iterator = 1; %tracks # of art2a iterations
    EndTest = 0; %stops script when == 1
    ClosestProximity = zeros(NumPID,1); 
    ClosestNeuron    = zeros(NumPID,1);
    LastNeuron = zeros(NumPID,1); 
    oldOrder = 1:NumPID;
    trackNeuron = []; 
    adjRate          = 1-LearningRate;  
  
%set up first particle  
    ClosestNeuron(1) = 1;
    NumClusters      = 1;
    RandomOrder = randperm(NumPID);  %get new ordering of 1:N
    RandomOrder = oldOrder(RandomOrder); %get new rand perm applied to old ordering
    [~, unperm] = sort(RandomOrder); %get 'undo' of rand perm
    oldOrder = RandomOrder; %save rand perm for next iteration
    WeightMatrix     = AreaMatrix(:,RandomOrder(1)); %for first iteration create a dummy weight matrix with first particle

%loop through algorithm
while EndTest == 0;
    
    %PFR track WM changes
    prevWM=WeightMatrix;
    
    %cluster particles
    for I = 1:NumPID 
        if mod(I,10000) == 0
            I
        end
            [ClosestProximity(I), ClosestNeuron(I)] = max(AreaMatrix(:,RandomOrder(I))'*WeightMatrix);
            if ClosestProximity(I) >= VigilanceFactor; %if particle matches a current cluster
                    temp = WeightMatrix(:,ClosestNeuron(I))*adjRate + (AreaMatrix(:,RandomOrder(I)) * LearningRate); %add to cluster and alter WM
                    WeightMatrix(:,ClosestNeuron(I)) = temp/norm(temp);
            else %if particle does not match any clusters
                NumClusters = NumClusters+1; 
                ClosestNeuron(I) = NumClusters;
                WeightMatrix = [WeightMatrix AreaMatrix(:,RandomOrder(I))]; %create new cluster
            end
    end
    
    %track WM
    if (size(WeightMatrix,2)==size(prevWM,2))  %ie no new clusters added
        ww = sum(WeightMatrix .* prevWM); %dot prod
        if (min(ww)>.9999)  %if all Wt vectors super close to previous then stop iterations early
            fprintf('INFO, Run Art2a, iteration: %i Weights converged \n',Iterator);
            EndTest=1;
        end;
    end;
    
    %find all PIDs that belong to each cluster after art2a iteration
    ClusterIDX = 1:NumClusters;
    UnpermNeuron = ClosestNeuron(unperm); %undo random perm
    
    %get 'tightness' metric for each cluster
    if avgP
        TempMAT = arrayfun(@(X) AreaMatrix(:,UnpermNeuron == X), ClusterIDX, 'UniformOutput',false);  %get area matrices for all particles in each cluster
        avgPROX{Iterator} = arrayfun(@(X) mean(TempMAT{X}'*WeightMatrix(:,X)), ClusterIDX); %get avg proximity of each particle to cluster WM
    end
    
    %find fraction of particles put into different clusters from last iteration
    if trackN
        ChangeNeuron = sum(UnpermNeuron ~= LastNeuron)/NumPID;
        trackNeuron = [trackNeuron ChangeNeuron]; %add fraction for this iteration to tracking metric
        LastNeuron = UnpermNeuron; %save results from this iteration for comparison later
        oldNumClusters = NumClusters;
    end
    
    %calculate PDistALL (note very time consuming)
    if pD
        D = arrayfun(@(X) pdist(AreaMatrix(:,UnpermNeuron == X)'), ClusterIDX, 'UniformOutput', false); %get average pair wise distance for all particles in a cluster
        Z = cellfun(@(X) squareform(X), D, 'UniformOutput',false);
        pDISTtemp = cellfun(@(X) mean(X(:)), Z);
        pDIST{Iterator} = pDISTtemp;
    end
    
    Iterator = Iterator + 1; %if reached MaxIteration end
    if Iterator > MaxIteration;
        EndTest = 1;
    end
    
    mem_stats=memory;
    fprintf('INFO, art2a, iter:%i,  tot mem used: %f Mb\n',Iterator,mem_stats.MemUsedMATLAB/1000000);
    
    %get new order of PIDs
    RandomOrder = randperm(NumPID);  %get new ordering of 1:N
    RandomOrder = oldOrder(RandomOrder); %get new rand perm applied to old ordering
    [~, unperm] = sort(RandomOrder); %get 'undo' of rand perm
    oldOrder = RandomOrder; %save rand perm for next iteration
    
end

%%
%sort the art2a results 

%get sorted PIDs in each cluster
outPartID = arrayfun(@(X) sortrows(inPartID(UnpermNeuron == X,:)), ClusterIDX, 'UniformOutput', false);
 
%remove empty cells and WM
IDnotEmpty = ~cellfun('isempty',outPartID);
outPartID = outPartID(IDnotEmpty);
WeightMatrix = WeightMatrix(:, IDnotEmpty);

%get clusters ordered from most to fewest # of particles
TempLength             = cellfun(@(X) size(X,1), outPartID, 'uniformoutput', false);
[partIDCount, SORTid] = sort(cell2mat(TempLength),'descend'); 
outPartID = outPartID(SORTid); %reorder clusters
outWeightMatrix = WeightMatrix(:,SORTid); %reorder weight matrix

% Truncate to get rid of all clusters that only contain 1 particle
idx = find(partIDCount > 1);
outPartID       = outPartID(idx);
partIDCount     = partIDCount(idx);
outWeightMatrix = outWeightMatrix(:,idx);

%put tracking metrics into matrices
if avgP
    avgPROXALL = padcat(avgPROX{1:min(Iterator,MaxIteration)}); %PFR added the iteration limit
    avgPROXALL = avgPROXALL(:,SORTid);
    avgPROXALL = avgPROXALL(:,idx);
end

if pD
    PDistALL = padcat(pDIST{:}); 
    PDistALL = PDistALL(:,SORTid); 
    PDistALL = PDistALL(:,idx); 
end



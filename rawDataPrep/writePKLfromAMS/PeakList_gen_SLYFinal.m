%% Method for baseline removal and peak finding of .ams spectra for ATOFMS.
%% Utilizes spline curves to create smooth lines, peak finding looks for
%% iterative max values. For more information see presentations and
%% descriptions on the ftp.
%% Jack Cahill - 3/14/13  -- version 1.4

% v1.2: Stupid LVN issues - instituted higher res baseline for small mz,
% refined for negative noise spikes and const for first 5 mz
% v1.3: BlownScale really changed things, redid how peaks are filtered and
% when areas were calculated. Saving are cals till the end
% v1.4: Instigated the while loop for peak duplicates. kept mz scan at 0.25
% but put const peak finding range of 1 mz.
% v2.0: Just added the ability to change the SN threshold
% v3.0: Re-did area calc

% input:
% % inData = full spectrum (just peaks) from .ams
% % IonType = positive or negative, positive = 0, negative = 1
% % calibIntercept/calibSlope = Calibration parameters (a = slope*(i-intercept^2))


% output:
%             FinalPeakSave(k,1) =  Height
%             FinalPeakSave(k,2) =  MZ
%             FinalPeakSave(k,3) =  Rel Area
%             FinalPeakSave(k,4) =  Final Area
%             FinalPeakSave(k,5) =  For area graph, min data pointer
%             FinalPeakSave(k,6) =  For area graph, max data pointer
%             FinalPeakSave(k,7) =  For area graph, peak data pointer
%             FinalPeakSave(k,8) =  Blow scale indicator (true if
%             height > 7500)

function [FinalPeakSave, Ring, inData, Noise, baseline, MZ] = PeakList_gen_SLYFinal(inData, IonType, calibIntercept, calibSlope, NoiseThresh)
% Call as [FinalPeakSave, Ring, inData, NoiseStd, baseline, MZ] = PeakList_gen_SLYFinal(inData, IonType, calibIntercept, calibSlope, NoiseThresh)
% Takes raw ams data (in time space) and calibration parameters and outputs calibrated peak picked data in m/z space.
% INPUTS
% inData: raw data output from get_spectrumAMS
% IonType: pos (0) or neg (1)
% calibIntercept: the intercept for the spectra calibration
% calibSlope - the slope for the spectra calibration
% NoiseThresh - S/N threshold to consider data, usually use 3 which is
% empirically derived   


NumPoints = length(inData); % should be 15000 but in case something changes

%%%%%%%%%%%%%%%%%% Baseline Calc %%%%%%%%%%%%%%%%%%%%%%

% Find local minima/median for every 10 mz range
% Fit a spline curve to these points (smooth polynomial)

%set up variables 
Start = calibSlope*(1-calibIntercept)^2; % finds the first mz bias
End = (calibSlope*(NumPoints-calibIntercept)^2); %last mz
MZ = (calibSlope*([1:NumPoints]-calibIntercept).^2)'; %all mz
Division = 25; % Mass division (not data point div!!) to find minima, make whole number, 10 was empirically derived to be the best.
j = 0; % serves as data point counter
mzdiv = Start:Division:200; %mz divisions

data_pointer1 = floor(sqrt(mzdiv/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
data_pointer2 = floor(sqrt((mzdiv+Division)/calibSlope)+calibIntercept);  % Finds ending MZ, in InData point space
data_pointer2(data_pointer2 > NumPoints) = NumPoints; %sanity checks, last data pointer can not exceed number of points
data_pointer1(data_pointer1 <= 0) = 1; %sanity checks, first data pointer cannot be <=0
localMin = zeros(length(mzdiv),2); % preallocate space for localMin, overestimate here

for i=Start:Division:200 % For every Division MZ: First 100 MZ needs to be higher resolution than latter MZ
     j=j+1;
    localMin(j,1) = min(inData(data_pointer1(j):data_pointer2(j))); % find min over those points
    localMin(j,2) = i+Division/2; % middle MZ position
end

% stops the begining from tailing wierdly
localMin = [localMin(1,1)*ones(4,1) [1:4]'; localMin];

% refine these values a bit
AvgBline(1,1) = mean(localMin(:,1));
AvgBline(1,2) = std(localMin(:,1));
inPeaktmp = localMin(:,1) > (AvgBline(1,1) + AvgBline(1,2)); % probably in a peak, so set to a slightly lower value
localMin(inPeaktmp,1) = AvgBline(1,1) + AvgBline(1,2);

%calc baseline for rest of spectra
Start = 200;
Division = 25; % Mass division (not data point div!!) to find minima, make whole number, 10 was empirically derived to be the best.
j = 0; % serves as data point counter
mzdiv = Start:Division:End;
data_pointer1 = floor(sqrt(mzdiv/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
data_pointer2 = floor(sqrt((mzdiv+Division)/calibSlope)+calibIntercept);  % Finds ending MZ, in InData point space
data_pointer2(data_pointer2 > NumPoints) = NumPoints; % Finds ending MZ, in InData point space
localMin2 = zeros(length(mzdiv),2);

for i=Start:Division:End % For every Division MZ
    j=j+1;
    localMin2(j,1) = min(inData(data_pointer1(j):data_pointer2(j))); % find min over those points
    localMin2(j,2) = i+Division/2; % middle MZ position
end

localMin = [localMin; localMin2]; %combine baseline calculated from beginning and end of spectra
baseline = spline(localMin(:,2),localMin(:,1),MZ); % fits a smooth curve to min points, output baseline of same size as data
tempPointer = floor(sqrt((5)/calibSlope)+calibIntercept); % Finds ending MZ, in InData point space
baseline(1:tempPointer) = localMin(1,1);
% plot(MZ,inData,'k',MZ,baseline,'r')

baseline(baseline < 0) = 0;  % can't have a baseline less than 0
inData = inData - baseline; % corrected data, overwrites previous data (saves memory)
inData(inData < 0) = 0; % can't have a negative value (just in case)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Noise Calculation  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Start = calibSlope*(1-calibIntercept)^2; % finds the first mz bias
Division = 5; % Mass division to find noise, Want larger than before so that it is less sensitive to peaks
j = 0; % serves as data point counter
mzdiv = Start:Division:200;
data_pointer1 = floor(sqrt(mzdiv/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
data_pointer2 = floor(sqrt((mzdiv+Division)/calibSlope)+calibIntercept);  % Finds ending MZ, in InData point space
data_pointer2(data_pointer2 > NumPoints) = NumPoints; %sanity checks, last data pointer can not exceed number of points
data_pointer1(data_pointer1 <= 0) = 1;  %sanity checks, first data pointer cannot be <=0
NoiseSave = zeros(length(mzdiv),3); % preallocate space for NoiseSave, overestimate here

for i=Start:Division:200 % For every Division mz
    j=j+1;
    inNoise = find(inData(data_pointer1(j):data_pointer2(j)) < 200 ); % Noise will never be higher than 200 after Bline correction, Find points < this value (empirical)
    inNoise = inNoise+data_pointer1(j)-1; % corrects pointer to match true data, side effect of using find
    NoiseSave(j,1) = std(inData(inNoise)); % Take std dev of those points    
    NoiseSave(j,3) = 0.5*median(inData(inNoise)); % Take median of those points
    % happens rarely but inNoise could be empty. Prob due to all data
    % points being above the 200 threshold
    if isnan(NoiseSave(j,1)) == 1
        if j ~= 1
            NoiseSave(j,1) = NoiseSave((j-1),1); % set to previous point, i.e. noise is static
        else
            NoiseSave(j,1) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        end
    end
    if isnan(NoiseSave(j,3)) == 1
        if j ~= 1
            NoiseSave(j,3) = NoiseSave((j-1),3); % set to previous point, i.e. median is static
        else
            NoiseSave(j,3) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        end
    end
    NoiseSave(j,2) = i+Division/2; % middle MZ position
end

Division = 10; % Mass division to find noise, Want larger than before so that it is less sensitive to peaks
mzdiv = 200:Division:End;
j = 0; % serves as data point counter
data_pointer1 = floor(sqrt(mzdiv/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
data_pointer2 = floor(sqrt((mzdiv+Division)/calibSlope)+calibIntercept);  % Finds ending MZ, in InData point space
data_pointer2(data_pointer2 > NumPoints) = NumPoints;
NoiseSave2 = zeros(length(mzdiv),3); % preallocate space for NoiseSave, overestimate here

for i=200:Division:End % For every Division mz
    j=j+1;
    inNoise = find(inData(data_pointer1(j):data_pointer2(j)) < 200 ); % Noise will never be higher than 200 after Bline correction, Find points < this value (empirical)
    inNoise = inNoise+data_pointer1(j)-1; % corrects pointer to match true data, side effect of using find
    NoiseSave2(j,1) = std(inData(inNoise)); % Take std dev of those points
    
    NoiseSave2(j,3) = 0.5*median(inData(inNoise)); % Take median of those points
    % happens rarely but inNoise could be empty. Prob due to all data
    % points being above the 200 threshold
    if isnan(NoiseSave2(j,1)) == 1
            NoiseSave2(j,1) = NoiseSave2((j-1),1); % set to previous point, i.e. noise is static
    end
    if isnan(NoiseSave2(j,3)) == 1
            NoiseSave2(j,3) = NoiseSave2((j-1),3); % set to previous point, i.e. median is static
    end    
    NoiseSave2(j,2) = i+Division/2; % middle MZ position
end

NoiseSave = [NoiseSave(1,1)*ones(4,1) [1:4]' zeros(4,1); NoiseSave; NoiseSave2]; %combine noise calculated from beginning and end of spectra, add beginning points to avoid weird tailing

for m=1:3  % # of iterations, totally empirical
    % % refine these values a bit
    AvgBline(1,1) = mean(NoiseSave(:,1));
    AvgBline(1,2) = std(NoiseSave(:,1));
    
    temp = find(NoiseSave(:,1) < (AvgBline(1,1)+AvgBline(1,2)));
    if isempty(temp) == 1
    else
        AvgBline(1,1) = mean(NoiseSave(temp,1));
        AvgBline(1,2) = std(NoiseSave(temp,1));
    end
    inPeaktmp = NoiseSave(:,1) > (AvgBline(1,1) + AvgBline(1,2)); % probably in a peak, so set to a slightly lower value
    NoiseSave(inPeaktmp,1) = AvgBline(1,1) + AvgBline(1,2);
end

NoiseStd = abs(spline(NoiseSave(:,2),NoiseSave(:,1),MZ)); % Combines median curve with std curve, abs() ensures no values < 0
NoiseMed = abs(spline(NoiseSave(:,2),NoiseSave(:,3),MZ)); % Combines median curve with std curve, abs() ensures no values < 0
NoiseStd = ceil(NoiseStd);

Noise(:,1) = NoiseStd;
Noise(:,2) = NoiseMed;

NoiseStd = NoiseThresh*NoiseStd + NoiseMed;
Noise(:,1) = NoiseStd;
% plot(MZ,inData,'k', MZ,Noise,'g')
NoiseThresh = 1; % reset to 1 because it has been incorporated into NoiseStd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PEAK FINDING %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%% PEAK FINDING %%%%%%%%%%%%%%%%%%%%%%%%
MZ_Division = 0.25; % Mass division to find peaks
j = 0; % serves as data point counter
BlowScale = 7500; % actuall around 8000 but after bline removal could make it be less
PeakSave = zeros(length(Start:MZ_Division:End),3);

%camille add
mzdiv = (Start:MZ_Division:(End-(3*MZ_Division)))'; %create mz divisions
data_pointer1 = (floor(sqrt((mzdiv/calibSlope))+calibIntercept));  %find lower location of mzDivision in data
data_pointer2 = (floor(sqrt(((mzdiv+MZ_Division)/calibSlope))+calibIntercept));  %find upper location of mzDivision in data
data_pointer1(data_pointer1 <= 0) = 1; %lower data location can't be <= 0
data_pointer2(data_pointer2 > NumPoints) = NumPoints; %upper data location can't ever be greater than NumPoints
lengthMZ = (data_pointer2-data_pointer1)+1; %find # of data points in each mz division
maxLength = max(lengthMZ); %find maximum # of data points

%divide data into a matrix to perform operations on
divideData = zeros(maxLength,length(data_pointer1));%variable to hold divided data
for i = 1:length(data_pointer1)
divideData(1:lengthMZ(i),i) = inData(data_pointer1(i):data_pointer2(i));
end
[localMax1 localMax2] = max(divideData,[],1); % find max over those points, #2 is index position
localMax1 = (localMax1');
localMax2 = (localMax2');
localMax2 = (localMax2+data_pointer1-1); % gets true position in inData
localMax2(localMax2 == NumPoints) = localMax2(localMax2 == NumPoints)-1; 

%find points above noise
sigNoise = localMax1./(NoiseStd(localMax2)); %find local max signal to noise ratio Note: Median is used here to estimate the true baseline hence std is added as well to noise in the calculation 
aboveNoise = sigNoise > NoiseThresh; %find points above S/N threshhold (using S/N thresh = 3 which is empirically derived
aboveNoiseIndex = find(aboveNoise); 
aboveNoiseLocalMax2 = localMax2(aboveNoise); %find position of points above SN thresh 
aboveNoiseLocalMax1 = localMax1(aboveNoise); %find value of points above SN thresh

% check that points to either side are real too so that blips arn't counted
if aboveNoiseLocalMax2(1)>1 
SN_BeforeNew = inData(aboveNoiseLocalMax2-1)./(NoiseStd(aboveNoiseLocalMax2-1));  
else
SN_BeforeNew = [0; inData(aboveNoiseLocalMax2(2:end)-1)./(NoiseStd(aboveNoiseLocalMax2(2:end)-1))]; % this is just for if a peak occurs right at the start of the spectra (ringing), makes it not a peak
end
SN_AfterNew = (inData(aboveNoiseLocalMax2+1)./(NoiseStd(aboveNoiseLocalMax2+1)));

%find real peaks
RealPeak = SN_BeforeNew < NoiseThresh & SN_AfterNew < NoiseThresh; % check that points to either side are real too so that blips arn't counted, don't do for points > 500 (assumed to be real)
RealPeak = RealPeak & aboveNoiseLocalMax1 < 250;
isPeak = RealPeak ~= 1; 
RealPeakIndex = aboveNoiseIndex(isPeak); %find location of real peaks
RealPeakIndexBlow = localMax1(RealPeakIndex) > BlowScale; %check if any peaks blow scale
RealPeakIndexNoBlow = localMax1(RealPeakIndex) < BlowScale; %location of peaks that don't blow scale

%initialize variables
numPeaks = length(RealPeakIndex);
UpperMZNew = zeros(numPeaks,1);
LowerMZNew = zeros(numPeaks,1);
upperPointerNew = zeros(numPeaks,1);
lowerPointerNew = zeros(numPeaks,1);
PeakSave = zeros(numPeaks,3);

%get data with peaks only (only need to deal with data considered peaks
%from now on
localMax1peak = localMax1(RealPeakIndex);
localMax2peak = localMax2(RealPeakIndex);
divideDataPeak = divideData(:,RealPeakIndex);
data_pointer1peak = data_pointer1(RealPeakIndex);
data_pointer2peak = data_pointer2(RealPeakIndex);

%treatment for blown peaks
if any(RealPeakIndexBlow)
    blowTempData = divideDataPeak(:,RealPeakIndexBlow); %find mz division where blown peaks are
    [blowTempRow, blowTempCol] = find(blowTempData > BlowScale); %get location of points that blow scale
    numBlow = size(blowTempData,2);
    RowMean = zeros(numBlow,1);
    RowMax = zeros(numBlow,1);
    RowMin = zeros(numBlow,1);
    for i = 1:numBlow % finds if there are multiple points that are blowing scale
        getRowData = blowTempRow(blowTempCol == i); %position of all points blowing scale
        RowMean(i,1) = mean(getRowData);
        RowMax(i,1) = max(getRowData);
        RowMin(i,1) = min(getRowData);
    end
    PeakSave(RealPeakIndexBlow,3) =  floor(RowMean)-1+data_pointer1peak(RealPeakIndexBlow);  % pointer position, roughly middle of blown peak
    % prep peak edge finding, puts lower and upper MZ at the
    % edges of the blown scale peak
    upperPointerNew(RealPeakIndexBlow) = RowMax-1+data_pointer1peak(RealPeakIndexBlow); % pointer in inData space
    lowerPointerNew(RealPeakIndexBlow) = RowMin-1+data_pointer1peak(RealPeakIndexBlow); % pointer in inData space  
    UpperMZNew(RealPeakIndexBlow) = (calibSlope*(upperPointerNew(RealPeakIndexBlow)-calibIntercept).^2) + MZ_Division; % look over the next division from max peak
    LowerMZNew(RealPeakIndexBlow) = (calibSlope*(lowerPointerNew(RealPeakIndexBlow)-calibIntercept).^2) - MZ_Division; % look over the next division from max peak
end

%no blown peaks
PeakSave(RealPeakIndexNoBlow,3) = localMax2peak(RealPeakIndexNoBlow); %save location of peak
UpperMZNew(RealPeakIndexNoBlow) = roundn(MZ(localMax2peak(RealPeakIndexNoBlow)),-2)+MZ_Division; % look over the next mz division from max peak
LowerMZNew(RealPeakIndexNoBlow) = roundn(MZ(localMax2peak(RealPeakIndexNoBlow)),-2)-MZ_Division; % look over the next mz division from max peak
LowerMZNew(LowerMZNew <= 0) = Start; %safety check, lowest mz can't be <= 0
lowerPointerNew(RealPeakIndexNoBlow) = localMax2peak(RealPeakIndexNoBlow)-1;
upperPointerNew(RealPeakIndexNoBlow) = localMax2peak(RealPeakIndexNoBlow)+1;

UpperMZNew = floor(sqrt(UpperMZNew/calibSlope)+calibIntercept); % convert to pointer
LowerMZNew = floor(sqrt(LowerMZNew/calibSlope)+calibIntercept); % convert to pointer
LowerMZNew(LowerMZNew<1) = 1; %safety checks
UpperMZNew(UpperMZNew>NumPoints) = NumPoints;

smallMZdiv = upperPointerNew > UpperMZNew | lowerPointerNew < LowerMZNew; %if the # of points per MZ_division is <3
if any(smallMZdiv)
    PeakSave(smallMZdiv,:) = [localMax2peak(smallMZdiv) localMax2peak(smallMZdiv) localMax2peak(smallMZdiv)]; %peak is a single point
end
normMZdiv = ~smallMZdiv; %if the # of points per MZ_division is >3 (should be almost all divisions)

%for all MZ divisions with > 3 points
normMZdivIndex = find(normMZdiv); %get index of all mz divisions with >3 points

%search over m/z division next to peaks to see if there is a greater value right next to it 
%set up variables
tempUpperPointer = upperPointerNew(normMZdiv);
tempUpperMZ = UpperMZNew(normMZdiv);
tempLowerPointer = lowerPointerNew(normMZdiv);
tempLowerMZ = LowerMZNew(normMZdiv);
lengthMZup = tempUpperMZ-tempUpperPointer+1;
lengthMZlow = tempLowerPointer-tempLowerMZ+1;
maxLengthUP = max(lengthMZup);
maxLengthLOW = max(lengthMZlow);
dataUP = nan(maxLengthUP,length(lengthMZup));
dataLOW = nan(maxLengthLOW,length(lengthMZlow));
indexUP = zeros(maxLengthUP,length(lengthMZup));
indexLOW = zeros(maxLengthLOW,length(lengthMZlow));
MaxID = zeros(length(normMZdivIndex),1);

for i = 1:length(normMZdivIndex)
    indexUPtemp = tempUpperPointer(i):tempUpperMZ(i); %index of data from point above peak to +mz division
    indexLOWtemp = tempLowerPointer(i):-1:tempLowerMZ(i); %index of data from point below peak to -mz division
    dataUP(1:lengthMZup(i),i) = inData(indexUPtemp); %data from point above peak to +mz division
    dataLOW(1:lengthMZlow(i),i) = inData(indexLOWtemp); %data from point below peak to -mz division
    indexUP(1:lengthMZup(i),i) = indexUPtemp; 
    indexLOW(1:lengthMZlow(i),i) =  indexLOWtemp;
end

belowThreshUP2 = bsxfun(@le, dataUP, 0.1*localMax1peak(normMZdivIndex)'); %find data less than threshold greater than +mz division
belowThreshUP3 = bsxfun(@le, dataUP,NoiseStd(localMax2peak(normMZdivIndex))');
belowThreshUP = double(belowThreshUP2 | belowThreshUP3);
belowThreshUP(isnan(dataUP)) = nan; 
belowThreshLOW2 = bsxfun(@le, dataLOW, 0.1*localMax1peak(normMZdivIndex)'); %find data less than threshold greater than -mz division
belowThreshLOW3 = bsxfun(@le, dataLOW,NoiseStd(localMax2peak(normMZdivIndex))');
belowThreshLOW = double(belowThreshLOW2 | belowThreshLOW3);
belowThreshLOW(isnan(dataLOW)) = nan;
newMaxUP = bsxfun(@gt, dataUP, localMax1peak(normMZdivIndex)'); %find data greater than current max for +mz division
newMaxLOW = bsxfun(@gt, dataLOW, localMax1peak(normMZdivIndex)'); %find data greater than current max for -mz division

belowThreshIDup = any(belowThreshUP,1); %find if any +mz division are below the threshold
belowThreshIDlow = any(belowThreshLOW,1); %find if any -mz division are below the threshold
newMaxIDUP = any(newMaxUP,1); %find if any +mz division are above current max
newMaxIDLOW = any(newMaxLOW,1); %find if any -mz division are above current max

converterUP = 0:maxLengthUP:maxLengthUP*(length(normMZdivIndex)-1); %row to index for matrix conversion
converterLOW = 0:maxLengthLOW:maxLengthLOW*(length(normMZdivIndex)-1); %row to index for matrix conversion

%for +m/z div: right edge of peak
if any(belowThreshIDup)
    rowUP = findfirst(belowThreshUP(:,belowThreshIDup)); %find the first index (row) where the value is below the threshhold
    edgeUP = rowUP+converterUP(belowThreshIDup); %convert the row to index
    PeakSave(normMZdivIndex(belowThreshIDup),2) = indexUP(edgeUP)+1; %save the right edge of the peak
    if any(newMaxIDUP(:,belowThreshIDup)) %check to see if any points exceeded the current max before falling below the threshold
        comboCheck = newMaxIDUP & belowThreshIDup; %find location of points exceeding the current max
        singleCheck = newMaxIDUP(belowThreshIDup);
        indexCombo = find(comboCheck);
        indexsingle = find(singleCheck);
        for j = 1:length(indexCombo); %for each mz division where the is a point exceeding the current max
            useData = newMaxUP(1:rowUP(indexsingle(j)),indexCombo(j)); %check to see if point happens before the R edge of the peak
            if any(useData) %if so 
            rowMax = find(useData,1,'last'); %find location of new peak
            convertMax = converterUP(indexCombo(j))+rowMax;
            PeakSave(normMZdivIndex(indexCombo(j)),3) = indexUP(convertMax); %save location of new peak
            MaxID(normMZdivIndex(indexCombo(j)))=1;
            end
        end
    end
end

if any(~belowThreshIDup) %if no R peak edge found
    PeakSave(normMZdivIndex(~belowThreshIDup),2) = tempUpperMZ(~belowThreshIDup); %right peak edge saved to edge of search (uppermz)
    if any(newMaxIDUP(:,~belowThreshIDup)) %check to see if any points exceeded the current max 
        comboCheck = newMaxIDUP & ~belowThreshIDup;  %find location of points exceeding the current max
        rowMax = findfirst(newMaxUP(:,comboCheck),1,1,'last');
        convertMax = converterUP(comboCheck)+rowMax;
        PeakSave(normMZdivIndex(comboCheck),3) = indexUP(convertMax); %save location of new peak
        MaxID(normMZdivIndex(comboCheck))=1;
    end
end

%for -m/z div: Left edge of peak
if any(belowThreshIDlow) %find the first index (row) where the value is below the threshhold
    rowLOW = findfirst(belowThreshLOW(:,belowThreshIDlow)); %find the first index (row) where the value is below the threshhold
    edgeLOW = rowLOW+converterLOW(belowThreshIDlow); %convert the row to index
    PeakSave(normMZdivIndex(belowThreshIDlow),1) = indexLOW(edgeLOW)-1;  %save the left edge of the peak
    if any(newMaxIDLOW(:,belowThreshIDlow)) %check to see if any points exceeded the current max before falling below the threshold
        comboCheck = newMaxIDLOW & belowThreshIDlow;  %find location of points exceeding the current max
        singleCheck = newMaxIDLOW(belowThreshIDlow);
        indexCombo = find(comboCheck);
        indexsingle = find(singleCheck);
        for j = 1:length(indexCombo);  %for each mz division where the is a point exceeding the current max
            useData = newMaxLOW(1:rowLOW(indexsingle(j)),indexCombo(j)); %check to see if point happens before the L edge of the peak
            if any(useData)
            rowMax = find(useData,1,'last');  %find location of new peak
            convertMax = converterLOW(indexCombo(j))+rowMax;
            PeakSave(normMZdivIndex(indexCombo(j)),3) = indexLOW(convertMax); %save location of new peak
            MaxID(normMZdivIndex(indexCombo(j)))=1;
            end
        end
    end
end

if any(~belowThreshIDlow)  %if no L peak edge found
    PeakSave(normMZdivIndex(~belowThreshIDlow),1) = tempLowerMZ(~belowThreshIDlow); %left peak edge saved to edge of search (uppermz)
    if any(newMaxIDLOW(:,~belowThreshIDlow))  %check to see if any points exceeded the current max 
        comboCheck = newMaxIDLOW & ~belowThreshIDlow; %find location of points exceeding the current max
        rowMax = findfirst(newMaxLOW(:,comboCheck),1,1,'last');
        convertMax = converterLOW(comboCheck)+rowMax;
        PeakSave(normMZdivIndex(comboCheck),3) = indexLOW(convertMax);  %save location of new peak
        MaxID(normMZdivIndex(comboCheck))=1;
    end    
end


% Saftey Check
if ~isempty(PeakSave) == 1
    PeakSave(PeakSave == 0) = 1;
    PeakSave(PeakSave > NumPoints) = NumPoints;
end

%for divisions where a new max was found
indexMaxID = find(MaxID);
if ~isempty(indexMaxID);
    for i=1:length(indexMaxID)
        lowID = PeakSave(indexMaxID(i),3):-1:PeakSave(indexMaxID(i),1);
        upID = PeakSave(indexMaxID(i),3):PeakSave(indexMaxID(i),2);
        newLedge = inData(lowID) <= 0.1*localMax1peak(indexMaxID(i)) | inData(lowID) <= NoiseStd(localMax2peak(indexMaxID(i))); %search for new left edge
        newRedge = inData(upID) <= 0.1*localMax1peak(indexMaxID(i)) | inData(upID) <= NoiseStd(localMax2peak(indexMaxID(i))); %search for new right edge
        if any(newLedge)
            PeakSave(indexMaxID(i),1) = lowID(find(newLedge,1))-1; %set new left edge
        end
        if any(newRedge)
            PeakSave(indexMaxID(i),2) = upID(find(newRedge,1))+1; %set new right edge
        end
    end
end

% Saftey Check
if ~isempty(PeakSave) == 1
    PeakSave(PeakSave == 0) = 1;
    PeakSave(PeakSave > NumPoints) = NumPoints;
end

%%
% Peak filtering - removing dublicates or merging overlaps
FinalPeakSave = [];
Ring = 0;
if isempty(PeakSave) == 1
    % This is if no peaks were found.  I have only found this to happen in
    % extremely ringy spectra - where the ring is bigger than any of the
    % actual peaks. In this case there is no way for any peaks to overcome
    % the noise threshold
    Ring = 1;
%         FinalPeakSave = [];
else
    
    % section of code will find duplicates, merge their pointers, and
    % truncate PeakSave to match true peak number
    % it helps to do this code first to greatly reduce the number of peaks
    Stop = length(PeakSave(:,1));
    %     for i = 1:length(moreRep)
    
    for i=1:Stop
        if Stop <= i
            break;
        end
        tempNum = find(PeakSave(:,1) <= PeakSave(i,3) & PeakSave(:,2) >= PeakSave(i,3));
        if length(tempNum) == 1 % found itself
        else % overlapping or matching peaks
            i=i-1;
            tempSave = zeros(length(PeakSave(:,1)),1);
            tempSave(tempNum(2:end)) = 1;
            PeakSave(tempNum(1),1) = min(PeakSave(tempNum,1));
            PeakSave(tempNum(1),2) = max(PeakSave(tempNum,2));
            PeakSave = PeakSave(~tempSave,:);
            Stop = length(PeakSave(:,1));
        end
        if Stop <= i
            break;
        end
    end
    
    % finds if Max pointer is same as min/max and checks that the peak
    % continues. Then it merges. It will iterate until, the number of peaks
    % stops changing
    k = 0;
    while k == 0
%         'inwhile'
        Stop = length(PeakSave(:,1));
        StopBefore = Stop; % num peaks used for comparision
        
        for i=1:Stop
            if Stop <= i
                break;
            end
            if PeakSave(i,1) == PeakSave(i,3) % max is on left edge
                tempPointer = PeakSave(i,1) - 1;
                if tempPointer > 0
                    if inData(tempPointer) > inData(PeakSave(i,1))
                        % Peak continues onward prob not a new peak
                        if i > 1
                            tempSave = zeros(length(PeakSave(:,1)),1);
                            tempSave(i) = 1;
                            PeakSave((i-1),2) = PeakSave(i,2); % sets the new max to combine
                            PeakSave = PeakSave(~tempSave,:);
                            Stop = length(PeakSave(:,1)); % change length
                        end
                    end
                end
            elseif PeakSave(i,2) == PeakSave(i,3) % max is on right edge
                tempPointer = PeakSave(i,2)+1;
                if tempPointer < NumPoints
                    if inData(tempPointer) > inData(PeakSave(i,2))
                        % peak continues
                        if i < Stop
                            tempSave = zeros(length(PeakSave(:,1)),1);
                            tempSave(i) = 1;
                            PeakSave((i+1),1) = PeakSave(i,1); % sets the new max to combine
                            PeakSave = PeakSave(~tempSave,:);
                            Stop = length(PeakSave(:,1)); % change length
                        end
                    end
                end
            end
            if Stop <= i
                break;
            end
        end
        
        % section of code will find duplicates, merge their pointers, and
        % truncate PeakSave to match true peak number
        Stop = length(PeakSave(:,1));
        for i=1:Stop
            if Stop <= i
                break;
            end
            tempNum = find(PeakSave(:,1) <= PeakSave(i,3) & PeakSave(:,2) >= PeakSave(i,3));
            %         tempNum = find(PeakSave(i,1) <= PeakSave(:,2) & PeakSave(i,2) >= PeakSave(:,1));
            if length(tempNum) == 1 % found itself
            else % overlapping or matching peaks
                i=i-1;
                tempSave = zeros(length(PeakSave(:,1)),1);
                tempSave(tempNum(2:end)) = 1;
                PeakSave(tempNum(1),1) = min(PeakSave(tempNum,1));
                PeakSave(tempNum(1),2) = max(PeakSave(tempNum,2));
                PeakSave = PeakSave(~tempSave,:);
                Stop = length(PeakSave(:,1));
            end
            if Stop <= i
                break;
            end
        end
        
        if Stop == StopBefore % no change, end the loop
            k = 1;
        end
    end
    
 %% now that recombination is over get new peak position
    lengthMZ = PeakSave(:,2)-PeakSave(:,1)+1; %find # of data points in each mz division
    maxLength = max(lengthMZ); %find maximum # of data points
    divideData = zeros(maxLength,length(lengthMZ)); %variable to hold divided data
    %divide data into mz divisions
    for i = 1:length(lengthMZ)
        divideData(1:lengthMZ(i),i) = inData(PeakSave(i,1):PeakSave(i,2));
    end
    [localMax1 localMax2] = max(divideData,[],1); % find max over those points, #2 is index position
    localMax2 = localMax2' + PeakSave(:,1) - 1; % gets true position in inData
    localMax1 = localMax1';
    noBlow = localMax1 <= BlowScale; % gets true position in inData
    PeakSave(noBlow,3) = localMax2(noBlow);
    if any(~noBlow)
        indexBlow = find(~noBlow);
        for i = 1:length(indexBlow)
            temp = find(inData(PeakSave(indexBlow(i),1):PeakSave(indexBlow(i),2)) > BlowScale); % finds if there are multiple points that are blowing scale
            PeakSave(indexBlow(i),3) = floor(mean(temp))-1+PeakSave(indexBlow(i),1); % pointer position, roughly middle of blown peak
        end
    end
    
    %%
    %remove points less than 5 mz and too close to end of spectra
    PeakSave = PeakSave(MZ(PeakSave(:,3))>5,:);
    PeakSave = PeakSave(PeakSave(:,2) <= (NumPoints-100),:);
    tempList = unique(MZ(PeakSave(:,3)),'first'); % finds unique MZ values
    NumList = length(tempList);
    if NumList == length(PeakSave(:,2))
    else
        tempPeakSave = zeros(size(PeakSave));
        k = 0;
        for i=1:NumList;
            tempNum = find(MZ(PeakSave(:,3)) == tempList(i));
            if length(tempNum) == 1
                k = k + 1;
                tempPeakSave(k,:) = PeakSave(tempNum,:); % no duplicate
            else
                k = k + 1;
                tempPeakSave(k,1) = min(PeakSave(tempNum,1)); % For area graph, min data pointer for area
                tempPeakSave(k,2) = max(PeakSave(tempNum,2)); % For area graph, max data pointer for area
                tempPeakSave(k,3) = max(PeakSave(tempNum,3)); % For area graph, max data pointer of peak
            end
        end
        PeakSave = tempPeakSave(1:k,:);
    end
    
    %% save data
    if ~isempty(PeakSave)
        FinalPeakSave(:,1) = (inData(PeakSave(:,3))); %height
        FinalPeakSave(:,2) = MZ(PeakSave(:,3),1); %mz
        FinalPeakSave(:,5:7) = PeakSave(:,:); %pointers
        %check for blow scale
        FinalPeakSave(:,8) = FinalPeakSave(:,1) > BlowScale;
        %calculate area
        for i=1:length(PeakSave(:,1))
            index = (PeakSave(i,1):PeakSave(i,2))';
            MZbin = ((calibSlope*((index+1)-calibIntercept).^2) - (calibSlope*((index-1)-calibIntercept).^2))/2; %get width of mz bins
            FinalPeakSave(i,4) = sum(inData(index).*MZbin); %multiply height by mz width
        end
        FinalPeakSave(:,3) = FinalPeakSave(:,4)/sum(FinalPeakSave(:,4)); %relative area
        if IonType == 1  % if a negative spectrum convert MZ to negative values           
            FinalPeakSave(:,2) = -FinalPeakSave(:,2);
        end
    else
    end
end
% Note: in yaada, you can find the number of peaks in a particle
% by: [a b c] = intercept(PEAK(:,'peakid'),PID(1))
end


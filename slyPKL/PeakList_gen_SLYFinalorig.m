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
% %             FinalPeakSave(k,1) =  Height
% %             FinalPeakSave(k,2) =  MZ
% %             FinalPeakSave(k,3) =  Rel Area
% %             FinalPeakSave(k,4) =  Final Area
% %             FinalPeakSave(k,5) =  For area graph, min data pointer
% %             FinalPeakSave(k,6) =  For area graph, max data pointer

function [FinalPeakSave Ring inData Noise baseline MZ] = PeakList_gen_SLYFinal(inData, IonType, calibIntercept, calibSlope, NoiseThresh)


NumPoints = length(inData); % should be 15000 but in case something changes
% new methods
% Find local minima/median for every 10 mz range
% Fit a spline curve to these points (smooth polynomial)
Start = calibSlope*(1-calibIntercept)^2; % finds the first mz bias
End = (calibSlope*(NumPoints-calibIntercept)^2);
for i=1:NumPoints
    MZ(i,1) = calibSlope*(i-calibIntercept)^2;  % Calc MZ for reference to later
end

Division = 25; % Mass division (not data point div!!) to find minima, make whole number, 10 was empirically derived to be the best.
j = 0; % serves as data point counter

for i=Start:Division:200 % For every Division MZ: First 100 MZ needs to be higher resolution than latter MZ
    j=j+1;
    data_pointer(1) = floor(sqrt(i/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
    data_pointer(2) = floor(sqrt((i+Division)/calibSlope)+calibIntercept); % Finds ending MZ, in InData point space
    if data_pointer(2) > NumPoints  % can't be greater than NumPoints, for last bin will be slightly smaller than normal
        data_pointer(2) = NumPoints;
    end
    localMin(j,1) = min(inData(data_pointer(1):data_pointer(2))); % find min over those points
    
    if j == 1 % stops the begining from tailing wierdly
        localMin(1:5,1) = localMin(1,1);
        localMin(1,2) = 1;
        localMin(2,2) = 2;
        localMin(3,2) = 3;
        localMin(4,2) = 4;
        j=j+4;
    end

    localMin(j,2) = i+Division/2; % middle MZ position
end

% refine these values a bit
AvgBline(1,1) = mean(localMin(:,1));
AvgBline(1,2) = std(localMin(:,1));
for i=1:length(localMin(:,1))
    if localMin(i,1) > (AvgBline(1,1) + AvgBline(1,2)) % probably in a peak, so set to a slightly lower value
        localMin(i,1) = AvgBline(1,1) + AvgBline(1,2);
    end
end

Start = 200;
Division = 25; % Mass division (not data point div!!) to find minima, make whole number, 10 was empirically derived to be the best.

for i=Start:Division:End % For every Division MZ
    j=j+1;
    data_pointer(1) = floor(sqrt(i/calibSlope)+calibIntercept); % Finds begining MZ, in InData point space
    data_pointer(2) = floor(sqrt((i+Division)/calibSlope)+calibIntercept); % Finds ending MZ, in InData point space
    if data_pointer(2) > NumPoints  % can't be greater than NumPoints, for last bin will be slightly smaller than normal
        data_pointer(2) = NumPoints;
    end
    localMin(j,1) = min(inData(data_pointer(1):data_pointer(2))); % find min over those points
    localMin(j,2) = i+Division/2; % middle MZ position
end

baseline = spline(localMin(:,2),localMin(:,1),MZ); % fits a smooth curve to min points, output baseline of same size as data
tempPointer = floor(sqrt((5)/calibSlope)+calibIntercept); % Finds ending MZ, in InData point space
baseline(1:tempPointer) = localMin(1,1);
% plot(MZ,inData,'k',MZ,baseline,'r')


for i=1:NumPoints % can't have a baseline less than 0
    if baseline(i) < 0
        baseline(i) = 0;
    end
    inData(i,1) = inData(i)-baseline(i); % corrected data, overwrites previous data (saves memory)
    if inData(i) < 0
        inData(i,1) = 0; % can't have a negative value (just in case)
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Noise Calculation  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Start = calibSlope*(1-calibIntercept)^2; % finds the first mz bias
Division = 5; % Mass division to find noise, Want larger than before so that it is less sensitive to peaks
j = 0; % serves as data point counter

temp = Start:Division:End;
NoiseSave = zeros((length(temp)+4),3); % preallocate space for localMin, overestimate here

for i=Start:Division:200 % For every Division mz
    j=j+1;
    data_pointer(1) = floor(sqrt(i/calibSlope)+calibIntercept); % Finds begining MZ
    data_pointer(2) = floor(sqrt((i+Division)/calibSlope)+calibIntercept); % Finds ending MZ
    if data_pointer(2) > NumPoints  % can't be greater than NumPoints, for last bin will be slightly smaller than normal
        data_pointer(2) = NumPoints;
    end
    inNoise = find(inData(data_pointer(1):data_pointer(2)) < 200 ); % Noise will never be higher than 200 after Bline correction, Find points < this value (empirical)
    inNoise = inNoise+data_pointer(1)-1; % corrects pointer to match true data, side effect of using find
    NoiseSave(j,1) = std(inData(inNoise)); % Take std dev of those points
    
    NoiseSave(j,3) = 0.5*median(inData(inNoise)); % Take median of those points
    % happens rarely but inNoise could be empty. Prob due to all data
    % points being above the 200 threshold
    if isnan(NoiseSave(j,1)) == 1
        if j == 1
            NoiseSave(j,1) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        else
            NoiseSave(j,1) = NoiseSave((j-1),1); % set to previous point, i.e. noise is static
        end
    end
    if isnan(NoiseSave(j,3)) == 1
        if j == 1
            NoiseSave(j,3) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        else
            NoiseSave(j,3) = NoiseSave((j-1),3); % set to previous point, i.e. median is static
        end
    end
    
    if j == 1
        NoiseSave(1:5,1) = NoiseSave(1,1);
        NoiseSave(1,2) = 1;
        NoiseSave(2,2) = 2;
        NoiseSave(3,2) = 3;
        NoiseSave(4,2) = 4;
        j=j+4;
    end

    NoiseSave(j,2) = i+Division/2; % middle MZ position
end



Division = 10; % Mass division to find noise, Want larger than before so that it is less sensitive to peaks
for i=200:Division:End % For every Division mz
    j=j+1;
    data_pointer(1) = floor(sqrt(i/calibSlope)+calibIntercept); % Finds begining MZ
    data_pointer(2) = floor(sqrt((i+Division)/calibSlope)+calibIntercept); % Finds ending MZ
    if data_pointer(2) > NumPoints  % can't be greater than NumPoints, for last bin will be slightly smaller than normal
        data_pointer(2) = NumPoints;
    end
    inNoise = find(inData(data_pointer(1):data_pointer(2)) < 200 ); % Noise will never be higher than 200 after Bline correction, Find points < this value (empirical)
    inNoise = inNoise+data_pointer(1)-1; % corrects pointer to match true data, side effect of using find
    NoiseSave(j,1) = std(inData(inNoise)); % Take std dev of those points
    
    NoiseSave(j,3) = 0.5*median(inData(inNoise)); % Take median of those points
    % happens rarely but inNoise could be empty. Prob due to all data
    % points being above the 200 threshold
    if isnan(NoiseSave(j,1)) == 1
        if j == 1
            NoiseSave(j,1) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        else
            NoiseSave(j,1) = NoiseSave((j-1),1); % set to previous point, i.e. noise is static
        end
    end
    if isnan(NoiseSave(j,3)) == 1
        if j == 1
            NoiseSave(j,3) = 0; % if first point just say noise/median = 0 It will adjust fast afterwards
        else
            NoiseSave(j,3) = NoiseSave((j-1),3); % set to previous point, i.e. median is static
        end
    end
    
    if j == 1
        NoiseSave(1:3,1) = NoiseSave(1,1);
        NoiseSave(1,2) = 1;
        NoiseSave(2,2) = 2;
        NoiseSave(3,2) = 3;
        NoiseSave(4,2) = 4;
        j=j+4;
    end
    
    NoiseSave(j,2) = i+Division/2; % middle MZ position
end

%remove empty rows of NoiseSave
NoiseSave = NoiseSave(1:j,:);


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
      
    for i=1:length(NoiseSave(:,1))
        if NoiseSave(i,1) > (AvgBline(1,1) + AvgBline(1,2)) % probably in a peak, so set to a slightly lower value
            NoiseSave(i,1) = AvgBline(1,1) + AvgBline(1,2);
        end
    end
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

MZ_Division = 0.25; % Mass division to find peaks
BlowScale = 7500; % actuall around 8000 but after bline removal could make it be less
j = 0; % serves as data point counter
PeakSave = zeros(length(Start:MZ_Division:End),3);
% PeakSave = [];
for i=Start:MZ_Division:(End-(3*MZ_Division)) % For every 1 (Division) mz, don't look at very end -- very noisy based on background calcs
    data_pointer(1) = floor(sqrt(i/calibSlope)+calibIntercept); % Finds begining MZ
    data_pointer(2) = floor(sqrt((i+MZ_Division)/calibSlope)+calibIntercept); % Finds ending MZ
    
    if data_pointer(2) > NumPoints  % can't be greater than NumPoints
        data_pointer(2) = NumPoints;
    end
    if data_pointer(1) <= 0  % can't be 0 or less
        data_pointer(1) = 1; % set to 1st point
    end
    [localMax(1) localMax(2)] = max(inData(data_pointer(1):data_pointer(2))); % find max over those points, #2 is index position
    localMax(2) = localMax(2) + data_pointer(1) - 1; % gets true position in inData
    if localMax(2) == NumPoints
        localMax(2) = localMax(2) -1;
    end
    
    if (localMax(1)/NoiseStd(localMax(2))) > NoiseThresh % if S/N > 3 then is a peak (3 is empirically derived), Note: Median is used here to estimate the true baseline hence std is added as well to noise in the calculation
        if (localMax(2)-1) < 1 % this is just for if a peak occurs right at the start of the spectra (ringing), makes it not a peak
            SN_Before = 0;
        else
            SN_Before = inData(localMax(2)-1)/(NoiseStd(localMax(2)-1)); % check that points to either side are real too so that blips arn't counted
        end
        SN_After = inData(localMax(2)+1)/(NoiseStd(localMax(2)+1));
        
        RealPeak = SN_Before < NoiseThresh && SN_After < NoiseThresh;
%         RealPeak = 0;
        RealPeak = RealPeak && localMax(1) < 250;
        if RealPeak == 1 % check that points to either side are real too so that blips arn't counted, don't do for points > 500 (assumed to be real)
            % Not a peak
        else % Is a peak
            j = j + 1; % peak Counter
            if localMax(1) > BlowScale
                % peak has blown scale, treat slightly different
                tempData = inData(data_pointer(1):data_pointer(2));
                temp = find(tempData > BlowScale); % finds if there are multiple points that are blowing scale
                
                PeakSave(j,3) = floor(mean(temp))-1+data_pointer(1); % pointer position, roughly middle of blown peak
                
                % prep peak edge finding, puts lower and upper MZ at the
                % edges of the blown scale peak
                upperPointer = max(temp)-1+data_pointer(1); % pointer in inData space
                lowerPointer = min(temp)-1+data_pointer(1); % pointer in inData space
                
                UpperMZ = (calibSlope*(upperPointer-calibIntercept)^2) + MZ_Division; % look over the next division from max peak, in this case 1 m/z
                LowerMZ = (calibSlope*(lowerPointer-calibIntercept)^2) - MZ_Division; % look over the next division from max peak, in this case 1 m/z
                if LowerMZ <= 0 % Saftey check
                    LowerMZ = Start;
                end
                
            elseif localMax(1) <= BlowScale
                % not blown scale
                
                PeakSave(j,3) = localMax(2);
                
                % Refine peak and calc areas - scan peak characteristics over
                % division space
                
                UpperMZ = roundn(MZ(localMax(2)),-2) + MZ_Division; % look over the next division from max peak, in this case 1 m/z
                LowerMZ = roundn(MZ(localMax(2)),-2) - MZ_Division; % look over the next division from max peak, in this case 1 m/z
                if LowerMZ <= 0 % Saftey check
                    LowerMZ = Start;
                end
                lowerPointer = localMax(2)-1;
                upperPointer = localMax(2)+1;
            end
            
            
            UpperMZ = floor(sqrt(UpperMZ/calibSlope)+calibIntercept); % convert to pointer
            LowerMZ = floor(sqrt(LowerMZ/calibSlope)+calibIntercept); % convert to pointer
            if LowerMZ < 1 % Saftey check
                LowerMZ = 1;
            end
            if UpperMZ > NumPoints % Saftey check
                UpperMZ = NumPoints;
            end
            
            maxChange = 0; % trigger if maximum was changed
            if upperPointer > UpperMZ % this can only happen if the # of points per MZ_Division is <3
                PeakSave(j,2) = localMax(2);
                PeakSave(j,3) = localMax(2);
            else
                for k=upperPointer:UpperMZ
                    if inData(k) <= 0.1*localMax(1) || inData(k) <= NoiseStd(localMax(2)) % stop if below noise or 10% of peak
                        PeakSave(j,2) = k+1; % position where peak ends to right
                        break
                    else
                        if inData(k) > localMax(1)
                            maxChange = 1;
                            PeakSave(j,3) = k;
                        end
                        if k == UpperMZ
                            PeakSave(j,2) = k; % never found an end, sumed everything
                        end
                    end
                end
            end
                
            if lowerPointer < LowerMZ  %  this can only happen if the # of points per MZ_Division is <3
                PeakSave(j,1) = localMax(2);
                PeakSave(j,3) = localMax(2);
            else
                
                for k=lowerPointer:-1:LowerMZ % Note: localMax(2)-1 so that the peak value isn't counted (already added)
                    if inData(k) <= 0.1*localMax(1) || inData(k) <= NoiseStd(localMax(2))
                        PeakSave(j,1) = k-1; % position where peak ends to left
                        break; % no need to look further
                    else
                        if inData(k) > localMax(1)
                            maxChange = 1;
                            PeakSave(j,3) = k; % set new max pointer
                        end
                        if k == LowerMZ
                            PeakSave(j,1) = k; % never found an end, sumed everything
                        end
                    end
                end      
            end
            
            if maxChange == 1
                for k=PeakSave(j,3):-1:PeakSave(j,1)
                    if inData(k) <= 0.1*localMax(1) || inData(k) <= NoiseStd(localMax(2))
                        PeakSave(j,1) = k-1; % position where peak ends to right
                        break; % no need to look further
                    end
                end
                for k=PeakSave(j,3):PeakSave(j,2)
                    if inData(k) <= 0.1*localMax(1) || inData(k) <= NoiseStd(localMax(2)) % stop if below noise or 10% of peak
                        PeakSave(j,2) = k+1; % position where peak ends to left
                        break
                    end
                end
            end
            
        end
    else % no peak so don't care
    end
end
% Remove Extras
PeakSave = PeakSave(1:j,:);

% Saftey Check
if isempty(PeakSave) == 1
else
    for i=1:length(PeakSave(:,1))
        for j=1:3
            if PeakSave(i,j) == 0
                PeakSave(i,j) = 1;
            elseif PeakSave(i,j) > NumPoints
                PeakSave(i,j) = NumPoints;
            end
        end
    end
end


% Peak filtering - removing dublicates or merging overlaps
Ring = 0;
FinalPeakSave = [];
if isempty(PeakSave) == 1
    % This is if no peaks were found.  I have only found this to happen in
    % extremely ringy spectra - where the ring is bigger than any of the
    % actual peaks. In this case there is no way for any peaks to overcome
    % the noise threshold
    %     FinalPeakSave = [];
else
    % section of code will find duplicates, merge their pointers, and
    % truncate PeakSave to match true peak number
    % it helps to do this code first to greatly reduce the peaks initially
    Stop = length(PeakSave(:,1));
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
        
        for i=1:length(PeakSave(:,1))
            [localMax(1) localMax(2)] = max(inData(PeakSave(i,1):PeakSave(i,2))); % find max over those points, #2 is index position
            localMax(2) = localMax(2) + PeakSave(i,1) - 1; % gets true position in inData
            if localMax(1) > BlowScale
                % peak has blown scale, treat slightly different
                temp = find(inData(PeakSave(i,1):PeakSave(i,2)) > BlowScale); % finds if there are multiple points that are blowing scale
                PeakSave(i,3) = floor(mean(temp))-1+PeakSave(i,1); % pointer position, roughly middle of blown peak
                
            elseif localMax(1) <= BlowScale
                % not blown scale
                PeakSave(i,3) = localMax(2);
            end
        end
        if Stop == StopBefore % no change, end the loop
            k = 1;
        end
    end
    
    
    % Now calc areas etc.
    for i=1:length(PeakSave(:,1))
        tempPeakSave(i,1) = inData(PeakSave(i,3)); % Height
        tempPeakSave(i,2) = MZ(PeakSave(i,3),1); % MZ
        
        % Calc area by normalizing the bin widths
        j = 0;
        Area = 0;
        for k=PeakSave(i,1):1:PeakSave(i,2)
            j = j+1;
            if k == 1
                MZbin = (calibSlope*((k+1)-calibIntercept)^2);  % Calc MZ for reference to later
            elseif k == NumPoints
                MZbin = (calibSlope*((k-1)-calibIntercept)^2);  % Calc MZ for reference to later
            else
                MZbin = ((calibSlope*((k+1)-calibIntercept)^2) - (calibSlope*((k-1)-calibIntercept)^2))/2;  % Calc MZ for reference to later
            end
            Area = Area + inData(k)*MZbin; % height is normalized to MZ division for these points
        end
        
        tempPeakSave(i,4) = Area; % Area
        tempPeakSave(i,5:7) = PeakSave(i,:); % pointers
    end
    
    % create final peak save
    tempList = unique(tempPeakSave(:,2),'first'); % finds unique MZ values
    NumList = length(tempList);
    k = 0;
    for i=1:NumList;
        tempNum = find(tempPeakSave(:,2) == tempList(i));
        if tempList(i) > 5 && max(tempPeakSave(tempNum,6)) <= (NumPoints-100) % has to be greater than m/z 5, if last data pointer  == ~End, then also false
            if length(tempNum) > 1
                k = k + 1;
                FinalPeakSave(k,1) = max(tempPeakSave(tempNum,1)); % Height
                FinalPeakSave(k,2) = max(tempPeakSave(tempNum,2)); % MZ
                %             FinalPeakSave(k,3) = floor(sqrt(FinalPeakSave(k,2)/calibSlope)+calibIntercept); % Finds data pointer, Pointer doesn't always match directly because area may be skewed to one side or the other
                FinalPeakSave(k,4) = max(tempPeakSave(tempNum,4)); % Final Area
                FinalPeakSave(k,5) = min(tempPeakSave(tempNum,5)); % For area graph, min data pointer
                FinalPeakSave(k,6) = max(tempPeakSave(tempNum,6)); % For area graph, max data pointer
                FinalPeakSave(k,7) = max(tempPeakSave(tempNum,7)); % For area graph, max data pointer
            else
                k = k + 1;
                FinalPeakSave(k,:) = tempPeakSave(tempNum,:); % no duplicate
            end
        end
    end
    

    if isempty(FinalPeakSave) == 1
    else
        FinalPeakSave(:,3) = FinalPeakSave(:,4)./sum(FinalPeakSave(:,4)); % calcs the relative peak area
        if IonType == 1  % if a negative spectrum convert MZ to negative values
            FinalPeakSave(:,2) = -FinalPeakSave(:,2);
        end
    end
end

if isempty(FinalPeakSave) == 1
    Ring = 1;
end

end


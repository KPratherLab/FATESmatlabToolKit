function varargout = calibFATES(varargin)
% Created with Matlab 2015a Guide by Gavin Cornwell, UCSD Mar 28 2016
% GUI for displaying and calibrating .ams files recorded by the ATOFMS
% Can manually control calibration values, edit table of points, load and
% save calibrations and other features.

% Note if using raw data files that are not .ams or .amz files will have to
% add your own code into loadSpectra to find the files (line 207) and 
% read in the file (229) 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calibFATES_OpeningFcn, ...
                   'gui_OutputFcn',  @calibFATES_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%------------------------------------Initialization---------------------------------------

function calibFATES_OpeningFcn(hObject, eventdata, handles, varargin) % initialization
handles.output = hObject;
guidata(hObject, handles);
global posMZ negMZ negA negB posA posB NumPoints particleSaveString calibToggle % declare global variables
    fu = mfilename('fullpath');% find folder for spectrumViewer
    NumPoints = 15000; % constant for number of points in spectrum
    calibToggle = 0; % calibration toggle set to default
    
    posA = 2E-6; % set default values for calibration parameters
    posB = 0; % set default values for calibration parameters
    negA = 2E-6; % set default values for calibration parameters
    negB = 0; % set default values for calibration parameters
    
    posMZ = 1:NumPoints; % uncalibrated spectrum 
    negMZ = 1:NumPoints; % uncalibrated spectrum
    for i = 1:20
        particleSaveString{i} = {''}; % create a cell of empty strings to be overwritten during the recording of which calibration points came from which particle
    end
    loadSpectra % read spectra from user input
    getSpectraIdx(handles) % determine list of spectra to toggle from--depends upon positive and negative checkboxes
    plotSpectra(handles,1,NumPoints) % plot spectra from data loaded by user

function varargout = calibFATES_OutputFcn(hObject, eventdata, handles) % required for GUI--put output here
varargout{1} = handles.output;

%-----------------------------------Operational functions----------------------
function navigateLeft(handles) % navigate to previous spectra
global filePointer posMZ negMZ foldDir polarity spectraIdx % declare global variables
    if isempty(filePointer) % make sure filePointer variable is real
        g = get(gca,'Title'); % determine particle name and then find its position in folder
        for i = 1:length(foldDir)
            if strcmp(foldDir(i).name,g) == 1
                filePointer = i; % find position of selected spectra in folder
                break
            end
        end
    end
    currentIdx = find(filePointer == spectraIdx); % find location of filePointer in list of spectra to be 
    if isempty(currentIdx) % this loop is to find the next filePointer if you are on a positive spectrum currently and then select only neg spectra 
        while 1
            filePointer = filePointer - 1;
            currentIdx = find(filePointer == spectraIdx);
            if ~isempty(currentIdx) > 0
                break
            end
        end
    end
    h = get(gca); % get handle for figure
    tmp = h.XLim; % get current xlim to save for drawing the next plot
    if polarity{filePointer} == 0 % positive
        tmpMin = abs(posMZ-tmp(1)); % find closest value for min
        tmpMax = abs(posMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    else % negative
        tmpMin = abs(negMZ-tmp(1)); % find closest value for min
        tmpMax = abs(negMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    end

    if ~isempty(filePointer)
        if filePointer == spectraIdx(1) % if at beginning of folder, don't redraw plot
        else % else redraw plot
            filePointer = spectraIdx(currentIdx-1); % go to next spectra
            plotSpectra(handles,newXMin,newXMax); % draw new spectrum
        end
    end
    
function navigateRight(handles) % navigate to later spectra
    global filePointer posMZ negMZ foldDir polarity spectraIdx% declare global variables
        if length(filePointer == 0)
        g = get(gca,'Title');
        for i = 1:length(foldDir)
            if strcmp(foldDir(i).name,g) == 1
                filePointer = i; % find position of selected spectra in folder
                break
            end
        end
    end
    currentIdx = find(filePointer == spectraIdx); % find location of filePointer in list of spectra to be 
    if isempty(currentIdx) % this loop is to find the next filePointer if you are on a positive spectrum currently and then select only neg spectra 
        while 1
            filePointer = filePointer + 1;
            currentIdx = find(filePointer == spectraIdx);
            if ~isempty(currentIdx) > 0
                break
            end
        end
    end
    h = get(gca); % get handle for figure
    tmp = h.XLim; % get current xlim to save for drawing the next plot
    if polarity{filePointer} == 0 % positive
        tmpMin = abs(posMZ-tmp(1)); % find closest value for min
        tmpMax = abs(posMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    else % negative
        tmpMin = abs(negMZ-tmp(1)); % find closest value for min
        tmpMax = abs(negMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    end

    if length(filePointer == 1)
        if filePointer == spectraIdx(end) % if at end of folder, don't redraw plot
        else % else redraw plot
            filePointer = spectraIdx(currentIdx+1); % go to next spectra
            plotSpectra(handles,newXMin,newXMax); % plot new spectrum
        end
    end
    
function figure1_KeyPressFcn(hObject, eventdata, handles) % navigate left and right
switch eventdata.Key
    case 'leftarrow' % navigate left
        navigateLeft(handles);
	case 'rightarrow' % navigate right
        navigateRight(handles);
end

function calculateMZ (NumPoints) % calculate MZ variables for x-axis
    global posMZ negMZ posA posB negA negB
    if isempty(posA) % assign default values if not already there
        posA = 2E-6;
    end
    if isempty(posB)
        posB = 0;
    end
    if isempty(negA)
        negA = 2E-6;
    end
    if isempty(negB)
        negB = 0;
    end
    NumPoints = 15000; % declare constant
    posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
    negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg 
        
function getSpectraIdx(handles) % determine index of particles for flipping through
global polIdx polarity spectraIdx
    polIdx = cell2mat(polarity); % convert cell to matrix for easy sorting
    spectraIdx = []; % declare blank variable
    if get(handles.positiveSpectraCheckBox,'Value') == 1; % if positive spectra box is checked
        spectraIdx = sort([spectraIdx find(polIdx == 0)]);
    end
    if get(handles.negativeSpectraCheckBox,'Value') == 1; % if negative spectra box is checked
        spectraIdx = sort([spectraIdx find(polIdx == 1)]);
    end
    
function loadSpectra % load spectrum from user input
    global filePointer defaultPath Data partSpeed partTime foldDir polarity NumPoints particleString filename % declare global variables
    [filename, pathName,fileIDX] = uigetfile({'*.ams','Open AMS file';'*.amz','Open AMZ file'},'Open AMS file',defaultPath); % user selects folder to load
    if ~filename
        error('Not a valid file name. Select a file to open') % error handling 
    end
    if fileIDX == 1
        foldDir = dir([pathName '*.ams']); % select ams files to read in
    elseif fileIDX == 2
        findAMZ = dir([pathName '*.amz']); % if no ams are present, check for amz
        if ~isempty(findAMZ)
            nameAMZ = {findAMZ.name}; %get pkl names
            nameAMZ(ismember(nameAMZ,{'avg_pos.amz', 'avg_neg.amz'})) = [];
            fullfileAMZ = fullfile(pathName,nameAMZ);
            if iscell(fullfileAMZ)
                fullfileAMS = cell(1,length(fullfileAMZ));
                for i = 1:length(fullfileAMZ)
                    tmp = unzip(fullfileAMZ{i},pathName); %uncompress amz file
                    fullfileAMS{i} = tmp{1};
                end
            else
                fullfileAMS = unzip(fullfileAMZ,pathName);
            end
        end
        foldDir = dir([pathName '*.ams']); % select ams files to read in
    elseif fileIDX == 3
        %if not using ams or amz files as raw tof files input 
        %add code here to create foldDir variable, which is a dir structure of raw data files in the
        %current folder that you would want to look at with calibFates
    end
    filePointer = 1; % create filePointer for keeping track of spectra 
    defaultPath = pathName; % create variable to keep track of folder spectra come from
    
    for i = 1:length(foldDir)
        if strcmp(foldDir(i).name,filename) == 1
            filePointer = i; % find position of selected spectra in folder
            break
        end
    end

    for i = 1:length(foldDir)
        %get_spectrumAMS
        if fileIDX < 3
        [~,polarity{i},partSpeed{i},partTime{i},~,Data{i}] = get_spectrumAMS(fullfile(pathName,foldDir(i).name)); % read in spectra
        fullfile(pathName,foldDir(i).name);
        Data{i} = Data{i}(1:NumPoints); % cut Data down to number of points
        partTime{i} = partTime{i}+datenum(2000,0,0)-1/24; % correct particle time
        else
            %if not using ams or amz files as raw tof files input
            %add code here to read raw data file (similar to get_spectrumAMS)
            %and then add data to variables
            %polarity (mass spec polarity), partSpeed (particle speed)
            %partTime (particle time), and Data (raw tof data) as shown
            %above
        end
        particleString{i} = [pathName foldDir(i).name]; % find full file path of each particle for calibration tracking
    end
    
    if fileIDX == 2
        delete(fullfile(pathName,'*.ams')); %delete uncompressed ams file
    end
    
function plotSpectra(handles,XMin,XMax) % plot current spectrum
    global filePointer polarity posMZ negMZ posA posB negA negB Data foldDir partSpeed partTime
    if polarity{filePointer} == 0 % load parameters for positive particles
        plot(posMZ,Data{filePointer})
        xlim([posMZ(XMin) posMZ(XMax)]) % set x limit
        set(handles.textA, 'String', num2str(posA),'FontSize',10); % display A
        set(handles.textB, 'String', num2str(posB),'FontSize',10); % display B
    else % load parameters for negative particles
        plot(negMZ,Data{filePointer})
        xlim([negMZ(XMin) negMZ(XMax)]) % set x limit
        set(handles.textA, 'String', num2str(negA),'FontSize',10); % display A
        set(handles.textB, 'String', num2str(negB),'FontSize',10); % display B
    end
  
    title(foldDir(filePointer).name); % title
    set(handles.textSpeed, 'String', num2str(partSpeed{filePointer})) % display particle speed
    set(handles.textDateTime, 'String', datestr(partTime{filePointer})) % diplay particle time
    ylim([0 max(Data{filePointer}(XMin:XMax))*1.05]); % reset ylim

% function [HitParticleCounter,IonType,Speed,TimeStampData,LaserPower,Data] = get_spectrumAMS(MS_filename) % function to read in AMS files--from Jack Cahill
% % read in binary file
%     filestream = fopen(MS_filename);
%     Version = fread(filestream,1,'short');
%     NumPoints = fread(filestream,1,'short');
%     HitParticleCounter = fread(filestream,1,'int');
%     ScatDelay = fread(filestream,1,'short');
%     Speed = fread(filestream,1,'float');
%     Size = fread(filestream,1,'float');
%     IonType = fread(filestream,1,'short');
%     TimeStamp_type = fread(filestream,1,'short');
%     TimeStampData = fread(filestream,1,'double');  % good
%     TimeText = fread(filestream,20,'*char');
% 
%     LaserPower = fread(filestream,1,'float');
%     ParticleHit = fread(filestream,1,'short');
%     FileType = fread(filestream,1,'short');
% 
%     FileName = (fread(filestream,80,'*char'));
% 
% % I assume all of these outputs are selectable from Tasware, like selecting
% % labels, re-caling and such
%     Cruft = fread(filestream,107,'*char');
%     TotalIntegral = fread(filestream,1,'int');
%     Baseline = fread(filestream,1,'short');
%     Calibrated = fread(filestream,1,'short');
%     CalibSlope = fread(filestream,1,'double');
%     CalibIntercept = fread(filestream,1,'double');
%     CalibData = fread(filestream,32,'double');
%     Labels = fread(filestream,1,'short');
%     Label = fread(filestream,480,'char'); % this value is actually a structure, 20 bytes of char, 4 bytes (float). Am not going to seperate
%     Cruft2 = fread(filestream,158, '*char');
% 
%     Data = fread(filestream,15000,'short');  % DATA!!
%     fclose(filestream);
    
function [A,B] = massCal(X,Y,plotFlag) % perform mass calibration with data
% massCal performs a calibration given two arrays of (1)
% points from a ATOFMS spectrum and (2) the corresponding M/Z values for
% those points. Adapted from WDR_DLL script (T Rebotier) 3/3/2016 using the
% following: The model is Y = A.(X-B)^2 where Y is the M/Z and X is the time (point number)
%We therefore seek the minimum of Q = SUM_i( (Y_i-A.(X_i-B)^2)^2 )
%for simplicity in the following formulas I omit the index i but keep in mind 
%that all sums are over different X_i's and Y_i's so only expressions in A or B purely
%can be factored out of the sums; we note Ui = (Xi-B)^2:
%Q = SUM( (Y-A.(X-B)^2)^2 )							(1)
%dQ/dA = 2 SUM ( (AU-Y)U )							(2)
%dQ/dA = 0 gives us A = SUM(YU)/SUM(U^2)				(3)
%dQ/dB = -4.A.SUM( (X-B).(AU-Y) )					(4)
%dQ/dB = 0 gives us B = SUM(X(AU-Y))/SUM(AU-Y)		(5)
%This formula cannot be solved directly but can be numerically (it converges)
%we iterate (3) and (5) starting with A =1 and B=0, after a few runs we keep the pair
%(a long and painful first order development of B(t+1) shows it 
%converges to the optimal B in the neighborhood of that optimal B...
%
%Changed (TR 7/22) to newB = average ( X - root(Y/A))
%A = slope, B = intercept
%To find M/Z using these two parameters:
% MZ = (A*([1:numPoints]-B).^2)';

if nargin == 2
    plotFlag = 0;
end

initial_A = 2E-6; % Define constants and initial values
initial_B = 0;
iterations = 1000;
A = initial_A;
B = initial_B;

for i = 1:iterations % Loop to determine A and B
    clear temp_B
    SYU = 0;SU2 = 0;newB = 0;
    U = (X-B).^2;
    U2 = U.*U;
    SYU = sum(Y.*U);
    SU2 = sum(U.*U);
    temp_B = X-(Y./A).^(0.5);
    A = SYU/SU2; % Define new A
    newB = sum(temp_B)/length(X); % Find new B, change slowly
    B = (0.9*B) + (0.1*newB);

end

%-----------------------------------GUI objects-----------------------------------------
function clearButton_Callback(hObject, eventdata, handles) % clear calibration table
set(handles.uitable1, 'Data', zeros(20,2)); % set cal table to be zeros

function calibrateButton_Callback(hObject, eventdata, handles) % Calibration button--temporary calibration
global filePointer foldDir Data NumPoints % declare global variables
    calData = get(handles.uitable1, 'Data'); % get data currently in table
    [A,B] = massCal(calData(calData(:,1)>0,1),calData(calData(:,1)>0,2)); % determine A and B from data in table
    tempMZ = A*([1:NumPoints]-B).^2'; % calculate temp MZ
    hold on % 
    plot(tempMZ,Data{filePointer}) % plot spectra with new MZ
    title(foldDir(filePointer).name); % display particle name as title
    hold off
    set(handles.textA, 'String', num2str(A),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(B),'FontSize',10); % display new B

function finalCalibrationButton_Callback(hObject, eventdata, handles) % final calibration
    global particleSaveString % declare global variables
    calData = get(handles.uitable1, 'Data'); % get data currently in table
    [A,B] = massCal(calData(calData(:,1)>0,1),calData(calData(:,1)>0,2)); % determine A and B from data in table
    calibrationPolarity = get(handles.negativeCalibrationButton,'Value')
    if calibrationPolarity == 0 % determine polarity by radio button
        [filename,pathName,~] = uiputfile('*.pcal','save calibration data as .pcal'); % user selects folder to load
    else
        [filename,pathName,~] = uiputfile('*.ncal','save calibration data as .ncal'); % user selects folder to load
    end
    fid = fopen([pathName filename],'w'); % open file from user input
    outFormat = '%s \n'; % output format
    fprintf(fid,outFormat,lastSpectrumPath); % write Path
    outFormat = '%d \n'; % output format
    fprintf(fid,outFormat,A); % write cal value
    fprintf(fid,outFormat,B); % write cal value
    for i = 1:length(calData) % write timepoint, mz, and particle it came from
        outFormat = '%d, %d, ';
        fprintf(fid,outFormat,calData(i,1),calData(i,2));
        outFormat = '%s \n';
        fprintf(fid,outFormat,char(particleSaveString{i}));
    end
    fclose(fid);
    set(handles.textA, 'String', num2str(A),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(B),'FontSize',10); % display new B
    
function toggleCalibrationButton_Callback(hObject, eventdata, handles) % button to toggle calibration on/off--go from uncalibrated to calibrated
global calibToggle NumPoints posA posB negA negB posMZ negMZ filePointer polarity
    h = get(gca);
    tmp = h.XLim;
    if polarity{filePointer} == 0 % positive
        tmpMin = abs(posMZ-tmp(1)); % find closest value for min
        tmpMax = abs(posMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    else % negative
        tmpMin = abs(negMZ-tmp(1)); % find closest value for min
        tmpMax = abs(negMZ-tmp(2)); % find closest value for max
        [~, newXMin] = min(tmpMin); %index of closest value
        [~, newXMax] = min(tmpMax); %index of closest value
    end
    if calibToggle == 1 % if toggle is on, turn it off
        posMZ = 1:15000;
        negMZ = 1:15000;
        calibToggle = 0;
    else % if toggle is off, turn it on
        calculateMZ (NumPoints);
        calibToggle = 1;
    end
    plotSpectra(handles, newXMin, newXMax) % plot new spectrum

function toggleCalibrationButton_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end
    
function positiveSpectraCheckBox_Callback(hObject, eventdata, handles) % determine spectraIdx
    getSpectraIdx(handles)
    
function negativeSpectraCheckBox_Callback(hObject, eventdata, handles) % determine spectraIdx
    getSpectraIdx(handles)
    
function positiveCalibrationButton_Callback(hObject, eventdata, handles) % ensure only positive or negative calibration is set
    fu = get(handles.positiveCalibrationButton,'Value');
    if fu == 1 
        set(handles.negativeCalibrationButton,'Value',0)
    else
        set(handles.negativeCalibrationButton,'Value',1)
    end
    
function negativeCalibrationButton_Callback(hObject, eventdata, handles) % ensure only positive or negative calibration is set
    fu = get(handles.negativeCalibrationButton,'Value');
    if fu == 1 
        set(handles.positiveCalibrationButton,'Value',0)
    else
        set(handles.positiveCalibrationButton,'Value',1)
    end
    
function uitable1_CreateFcn(hObject, eventdata, handles) % create table for holding calibration values

function uitable1_CellEditCallback(hObject, eventdata, handles) % make sure particle data is recorded when editing calibration table
global filePointer particleString particleSaveString % declare global variables
    [X] = eventdata.Indices; % get new value put int
    particleSaveString{X(1)} = particleString{filePointer}; % get string for current particle and save it

function figure1_DeleteFcn(hObject, eventdata, handles) % clear global variables from workspace upon exiting
clear -global filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints...
    particleString particleSaveString exportPath % clear global variables

%--------------------------Menu and tooolbar button functions---------------------------
function fileMenu_Callback(hObject, eventdata, handles) % create menu->file

function viewMenu_Callback(hObject, eventdata, handles) % create menu->view

function calibrationMenu_Callback(hObject, eventdata, handles) % create menu->calibration    

function openSpectrum_Callback(hObject, eventdata, handles) % open new spectrum code, from File menu
global NumPoints % declare global variables
    loadSpectra % read spectra from user input
    plotSpectra(handles,1,NumPoints) % plot spectra from data loaded by user

function exitFigure_Callback(hObject, eventdata, handles) % quit program and clear global variables
clear -global filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints  particleString particleSaveString exportPath % clear global variables
    close(gcf) % close figure
    
function resetView_Callback(hObject, eventdata, handles) % reset view of plot normal
    global filePointer Data posMZ negMZ polarity NumPoints % declare global variables
    if polarity{filePointer} == 0 % determine if particle is positive or negative
        xlim([1 posMZ(NumPoints)]); % reset to posMZ upper limit
    else
        xlim([1 negMZ(NumPoints)]); % reset to negMZ upper limit
    end
    ylim([0 max(Data{filePointer})*1.05]); % reset y-axis

function setX_Callback(hObject, eventdata, handles) % set X axis limts manually
    temp_min = inputdlg('Enter minimum X-value'); % min value input
    temp_max = inputdlg('Enter maximum X-value'); % max value input 
    xlim([str2num(char(temp_min)) str2num(char(temp_max))]) % set limits from user inputs

function setY_Callback(hObject, eventdata, handles) % set Y axis limts manually
    temp_min = inputdlg('Enter minimum Y-value'); % min value input
    temp_max = inputdlg('Enter maximum Y-value'); % max value input
    ylim([str2num(char(temp_min)) str2num(char(temp_max))]) % set limits from user inputs

function setCalibrationValues_Callback(hObject, eventdata, handles) % set A and B calibration values manually
global filePointer posMZ negMZ negA negB posA posB polarity NumPoints calibToggle% declare global variables
    tempA = inputdlg('Enter new A value'); % user input for A
    tempB = inputdlg('Enter new B value'); % user input for B
    if polarity{filePointer} == 0 % positive
        posA = str2num(char(tempA)); % reset A
        posB = str2num(char(tempB)); % reset B
        posMZ = (posA*([1:NumPoints]-posB).^2)'; % recalculate posMZ array
        set(handles.textA, 'String', num2str(posA),'FontSize',10); % display new A
        set(handles.textB, 'String', num2str(posB),'FontSize',10); % display new B
        plotSpectra(handles,1,posMZ(NumPoints))
    else % negative
        negA = str2num(char(tempA)); % user input for A
        negB = str2num(char(tempB)); % user input for B
        negMZ = (negA*([1:NumPoints]-negB).^2)'; % recalculate negMZ array
        set(handles.textA, 'String', num2str(negA),'FontSize',10); % display new A
        set(handles.textB, 'String', num2str(negB),'FontSize',10); % display new B
        plotSpectra(handles,1,negMZ(NumPoints))
    end
    calibToggle = 1;

function addPoint_Callback(hObject, eventdata, handles) % add point for calibration--very important function
    global filePointer posMZ negMZ polarity particleString particleSaveString % declare global variables
    [X,height] = ginput(1); % activate user input through mouse
    if polarity{filePointer} == 0
        tmp = abs(posMZ-X);
        [~, idx] = min(tmp); %index of closest value
    else
        tmp = abs(negMZ-X);
        [~, idx] = min(tmp); %index of closest value
    end
    set(handles.textPoints, 'String', sprintf('%s, %s \n %s',num2str(X),num2str(idx),num2str(height))); % output index and value to textbox
    [Y] = inputdlg('Please enter M/Z'); % get user input MZ for point

    calData = get(handles.uitable1, 'Data'); % read current table
    calCounter = find(calData(:,1) == 0); % find first value with zero for the timepoint value
        if polarity{filePointer} == 0 % find particle polarity
            tmp = abs(posMZ-X); % find closest value
            [~, idx] = min(tmp); %index of closest value
            particleSaveString{calCounter(1)} = particleString{filePointer}; % save particle name for writing to the calibration file
            calData(calCounter(1),1) = idx; % fill in first empty value 
            calData(calCounter(1),2) = str2num(char(Y)); % fill in first empty value
            set(handles.uitable1, 'Data', [calData]); % write to table           
        else
            tmp = abs(negMZ-X); % find closest value
            [~, idx] = min(tmp); %index of closest value
            particleSaveString{calCounter(1)} = particleString{filePointer}; % save particle name for writing to the calibration file
            calData(calCounter(1),1) = idx; % fill in first empty value 
            calData(calCounter(1),2) = str2num(char(Y)); % fill in first empty value
            set(handles.uitable1, 'Data', [calData]); % write to table
        end
    set(handles.textPoints, 'String', sprintf('%s, %s',num2str(X),num2str(idx))); % output

function getPoint_Callback(hObject, eventdata, handles) % get point and display current mz and timepoint
    global polarity filePointer posMZ negMZ
    [X,height] = ginput(1); % user input prompt
    if polarity{filePointer} == 0 % positive
        tmp = abs(posMZ-X);
        [~, idx] = min(tmp); %index of closest value
    else % negative
        tmp = abs(negMZ-X);
        [~, idx] = min(tmp); %index of closest value
    end
    set(handles.textPoints, 'String', sprintf('%s, %s \n %s',num2str(X),num2str(idx),num2str(height))); % output index and value to textbox

function loadCalibration_Callback(hObject, eventdata, handles) % load user-saved calibration file
    global filePointer negA negB posA posB polarity
    calibrationPolarity = get(handles.negativeCalibrationButton,'Value');
    if polarity{filePointer}  == 0
        [filename,pathName,~] = uigetfile('*.pcal','load calibration data for positive spectra'); % user selects folder to load
    else
        [filename,pathName,~] = uigetfile('*.ncal','load calibration data for negative spectra'); % user selects folder to load
    end
    fid = fopen([pathName filename],'r');
    fgetl(fid);
    A = str2num(fgetl(fid)) % load cal value
    B = str2num(fgetl(fid)) % load cal value
    fclose(fid);
    polarity{filePointer}
    if polarity{filePointer}  == 0
        posA = A;
        posB = B;
    else
        negA = A;
        negB = B;
    end
    set(handles.textA, 'String', num2str(posA),'FontSize',10); % display A
    set(handles.textB, 'String', num2str(posB),'FontSize',10); % display B
    
function saveCalibration_Callback(hObject, eventdata, handles)% save CURRENT calibration values from the text boxes below the plot
global defaultPath polarity filePointer % clear global variables
    A = str2num(get(handles.textA, 'String')); % display new A
    B = str2num(get(handles.textB, 'String')); % display new B
    calibrationPolarity = get(handles.negativeCalibrationButton,'Value')
    if polarity{filePointer}  == 0
        [filename,pathName,~] = uiputfile('*.pcal','load calibration data for positive spectra'); % user selects folder to load
    else
        [filename,pathName,~] = uiputfile('*.ncal','load calibration data for negative spectra'); % user selects folder to load
    end

    fid = fopen([pathName filename],'w');
    outFormat = '%s \n'; % output format
    fprintf(fid,outFormat,[defaultPath]); % write Path
    outFormat = '%d \n'; % output format
    fprintf(fid,outFormat,A); % write cal value
    fprintf(fid,outFormat,B); % write cal value
    fclose(fid)

function leftArrowButton_ClickedCallback(hObject, eventdata, handles) % arrow button to navigate to previous spectrum
    navigateLeft(handles)
    
function rightArrowButton_ClickedCallback(hObject, eventdata, handles) % arrow button to navigate to next spectrum
    navigateRight(handles)

function clearButton_KeyPressFcn(hObject, eventdata, handles) % clear calibration table
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function calibrateButton_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function finalCalibrationButton_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function positiveSpectraCheckBox_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function negativeSpectraCheckBox_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function positiveCalibrationButton_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

function negativeCalibrationButton_KeyPressFcn(hObject, eventdata, handles) % still go left and right even if focus is on toggle calibration button
    switch eventdata.Key
        case 'leftarrow' % navigate left
            navigateLeft(handles);
        case 'rightarrow' % navigate right
            navigateRight(handles);
    end

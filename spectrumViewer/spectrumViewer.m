function varargout = spectrumViewer(varargin)
% Created with Matlab 2015a Guide by Gavin Cornwell, UCSD Mar 28 2016
% GUI for displaying and calibrating .ams files recorded by the ATOFMS
% Can manually control calibration values, edit table of points, load and
% save calibrations and other features.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spectrumViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @spectrumViewer_OutputFcn, ...
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

function spectrumViewer_OpeningFcn(hObject, eventdata, handles, varargin) % initialization
handles.output = hObject;
guidata(hObject, handles);

global filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints particleString particleSaveString exportPath% declare global variables
fu = mfilename('fullpath');% find folder for spectrumViewer
exportPath = fu(1:end-14);
fid = fopen([exportPath 'defaultConfig.cal']); % load config values
defaultPath = fgetl(fid); % load Path 
posA = str2num(fgetl(fid)); % load cal value
posB = str2num(fgetl(fid)); % load cal value
negA = str2num(fgetl(fid)); % load cal value
negB = str2num(fgetl(fid)); % load cal value
fclose(fid);

NumPoints = 15000; % declare constant
posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg
for i = 1:20
    particleSaveString{i} = {''}; % create a cell of empty strings to be overwritten during the recording of which calibration points came from which particle
end

[filename, pathName,~] = uigetfile({'*.ams','Open AMS file';'*.amz','Open AMZ file'},'Open AMS file',defaultPath); % user selects folder to load
if ~filename
    error('select an AMS file to open') % error handling so that not selecting a file does not overwrite your default configurations
end

foldDir = dir([pathName '*.ams']); % select ams files to read in
if isempty(foldDir )
    foldDir = dir([pathName '*.amz']);
end
filePointer = 1; % create filePointer for keeping track of spectra
defaultPath = pathName; % create variable to keep track of folder spectra come from

fu = mfilename('fullpath'); % find folder for spectrumViewer
fid = fopen([fu(1:end-20) 'defaultConfig.cal'],'w'); %load file for writing
[fu(1:end-20) 'defaultConfig.cal'];
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,posA); % write cal value
fprintf(fid,outFormat,posB); % write cal value
fprintf(fid,outFormat,negA); % write cal value
fprintf(fid,outFormat,negB); % write cal value
fclose(fid);

for i = 1:length(foldDir)
    if strcmp(foldDir(i).name,filename) == 1
        filePointer = i; % find position of selected spectra in folder
        break
    end
end

for i = 1:length(foldDir)
    [~,polarity{i},partSpeed{i},partTime{i},~,Data{i}] = get_spectrumAMS(fullfile(pathName,foldDir(i).name)); % read in spectra
    fullfile(pathName,foldDir(i).name);
    Data{i} = Data{i}(1:NumPoints); % cut Data down to number of points
    partTime{i} = partTime{i}+datenum(2000,0,0)-1/24; % correct particle time
    particleString{i} = [pathName foldDir(i).name]; % find full file path of each particle for calibration tracking
end

% plot figure
if polarity{filePointer} == 0 % load parameters for positive particles
    plot(posMZ,Data{filePointer})
    xlim([0 posMZ(NumPoints)]) % set x limit
    set(handles.textA, 'String', num2str(posA),'FontSize',10); % display A
    set(handles.textB, 'String', num2str(posB),'FontSize',10); % display B
else % load parameters for negative particles
    plot(negMZ,Data{filePointer})
    xlim([0 negMZ(NumPoints)]) % set x limit
    set(handles.textA, 'String', num2str(negA),'FontSize',10); % display A
    set(handles.textB, 'String', num2str(negB),'FontSize',10); % display B
end

title(foldDir(filePointer).name); % display particle name as title
set(handles.textSpeed, 'String', num2str(partSpeed{filePointer})) % display particle speed
ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % set ylim to nominal baseline and just above max Y
set(handles.textDateTime, 'String', datestr(partTime{filePointer})) % display particle

function varargout = spectrumViewer_OutputFcn(hObject, eventdata, handles) % required for GUI--put output herre
varargout{1} = handles.output;

%% Keyboard navigation
function figure1_KeyPressFcn(hObject, eventdata, handles) % function for cycling through spectra
global filePointer Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity% declare global variables
switch eventdata.Key
    case 'leftarrow' % navigate left
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
        if filePointer == 1 % if at beginning of folder, don't redraw plot
        else % else redraw plot
            filePointer = filePointer-1; % go to next spectra
            if polarity{filePointer} == 0 % positive
                plot(posMZ,Data{filePointer}) % redraw
                xlim([posMZ(newXMin) posMZ(newXMax)]) % set xlim as determined previously
                set(handles.textA, 'String', num2str(posA),'FontSize',10); % display correct cal values
                set(handles.textB, 'String', num2str(posB),'FontSize',10);
            else % negative
                plot(negMZ,Data{filePointer}) % go to next spectra
                xlim([negMZ(newXMin) negMZ(newXMax)]) % set xlim to keep timepoint the same
                set(handles.textA, 'String', num2str(negA),'FontSize',10); % display correct cal values
                set(handles.textB, 'String', num2str(negB),'FontSize',10);
            end
            title(foldDir(filePointer).name); % title
            set(handles.textSpeed, 'String', num2str(partSpeed{filePointer})) % display particle speed
            set(handles.textDateTime, 'String', datestr(partTime{filePointer})) % diplay particle time
            ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % reset ylim
        end
	case 'rightarrow' % navigate right
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
        if filePointer == length(foldDir)
        else % else redraw plot
            filePointer = filePointer+1; % go to next spectra
            if polarity{filePointer} == 0 % positive
                plot(posMZ,Data{filePointer}) % redraw
                xlim([posMZ(newXMin) posMZ(newXMax)]) % set xlim as determined previously
                set(handles.textA, 'String', num2str(posA),'FontSize',10); % display correct cal values
                set(handles.textB, 'String', num2str(posB),'FontSize',10);
            else % negative
                plot(negMZ,Data{filePointer}) % go to next spectra
                xlim([negMZ(newXMin) negMZ(newXMax)]) % set xlim to keep timepoint the same
                set(handles.textA, 'String', num2str(negA),'FontSize',10); % display correct cal values
                set(handles.textB, 'String', num2str(negB),'FontSize',10);
            end
            title(foldDir(filePointer).name); % title
            set(handles.textSpeed, 'String', num2str(partSpeed{filePointer})) % display particle speed
            set(handles.textDateTime, 'String', datestr(partTime{filePointer})) % diplay particle time
            ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % reset ylim
        end
end

 %% Button
function clearButton_Callback(hObject, eventdata, handles) % clear calibration table
set(handles.uitable1, 'Data', zeros(20,2)); % set cal table to be zeros

function calibrateButton_Callback(hObject, eventdata, handles)
global filePointer foldDir Data NumPoints % declare global variables
calData = get(handles.uitable1, 'Data'); % get data currently in table
[A,B] = massCalibration(calData(calData(:,1)>0,1),calData(calData(:,1)>0,2)); % determine A and B from data in table
tempMZ = A*([1:NumPoints]-B).^2'; % calculate temp MZ
hold on % 
plot(tempMZ,Data{filePointer}) % plot spectra with new MZ
title(foldDir(filePointer).name); % display particle name as title
hold off
set(handles.textA, 'String', num2str(A),'FontSize',10); % display new A
set(handles.textB, 'String', num2str(B),'FontSize',10); % display new B

function finalCalibrationButton_Callback(hObject, eventdata, handles)
global filePointer defaultPath exportPath partTime negA negB posA posB polarity particleSaveString % declare global variables
A = str2num(get(handles.textA,'String'));
B = str2num(get(handles.textB,'String'));
if polarity{filePointer} == 0
    posA = A;
    posB = B;
    polMarker = 'Positive';
else
    negA = A;
    negB = B;
    polMarker = 'Negative';
end

fu = mfilename('fullpath'); % find folder for spectrumViewer
fid = fopen([fu(1:end-20) 'defaultConfig.cal'],'w'); %load file for writing
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,posA); % write cal value
fprintf(fid,outFormat,posB); % write cal value
fprintf(fid,outFormat,negA); % write cal value
fprintf(fid,outFormat,negB); % write cal value
fclose(fid);

filepath = datestr(partTime{1},'yyyymmdd');
mkdir(sprintf('C:/calibration/%s',filepath));
filename = sprintf('C:/calibration/%s/cal_%s.cal',filepath,datestr(datetime('now'),'yyyymmddHHMMSS'));
calData = get(handles.uitable1, 'Data'); % get data currently in table
fid = fopen(filename,'w');
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,posA); % write cal value
fprintf(fid,outFormat,posB); % write cal value
fprintf(fid,outFormat,negA); % write cal value
fprintf(fid,outFormat,negB); % write cal value
for i = 1:length(calData)
    outFormat = '%d, %d, ';
    fprintf(fid,outFormat,calData(i,1),calData(i,2));
    outFormat = '%s \n';
    fprintf(fid,outFormat,char(particleSaveString{i}));
end
fclose(fid);

%% Menu creation and commands
function uitable1_CreateFcn(hObject, eventdata, handles) % create table for holding calibration values

function fileMenu_Callback(hObject, eventdata, handles) % create menu->file

function viewMenu_Callback(hObject, eventdata, handles) % create menu->view

function calibrationMenu_Callback(hObject, eventdata, handles) % create menu->calibration

function openSpectrum_Callback(hObject, eventdata, handles) % open new spectrum code, from File menu
global filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints particleString exportPath% declare global variables

NumPoints = 15000; % declare constant
posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg

NumPoints = 15000; % declare constant
posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg

[filename, pathName,~] = uigetfile({'*.ams','Open AMS file';'*.amz','Open AMZ file'},'Open AMS file',defaultPath); % user selects folder to load
if ~filename
    error('select an AMS file to open') % error handling so that not selecting a file does not overwrite your default configurations
end

foldDir = dir([pathName '*.ams']); % select ams files to read in
filePointer = 1; % create filePointer for keeping track of spectra
defaultPath = pathName; % create variable to keep track of folder spectra come from

fu = mfilename('fullpath'); % find folder for spectrumViewer
fid = fopen([fu(1:end-20) 'defaultConfig.cal'],'w'); %load file for writing
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,posA); % write cal value
fprintf(fid,outFormat,posB); % write cal value
fprintf(fid,outFormat,negA); % write cal value
fprintf(fid,outFormat,negB); % write cal value
fclose(fid);

for i = 1:length(foldDir)
    if strcmp(foldDir(i).name,filename) == 1
        filePointer = i; % find position of selected spectra in folder
        break
    end
end

for i = 1:length(foldDir)
    [~,polarity{i},partSpeed{i},partTime{i},~,Data{i}] = get_spectrumAMS(fullfile(pathName,foldDir(i).name)); % read in spectra
    Data{i} = Data{i}(1:NumPoints); % cut Data down to number of points
    partTime{i} = partTime{i}+datenum(2000,0,0)-1/24; % correct particle time
    particleString{i} = [pathName foldDir(i).name]; % find full file path of each particle for calibration tracking
end

% plot figure
if polarity{filePointer} == 0 % load parameters for positive particles
    plot(posMZ,Data{filePointer})
    xlim([0 posMZ(NumPoints)]) % set x limit
    set(handles.textA, 'String', num2str(posA),'FontSize',10); % display A
    set(handles.textB, 'String', num2str(posB),'FontSize',10); % display B
else % load parameters for negative particles
    plot(negMZ,Data{filePointer})
    xlim([0 negMZ(NumPoints)]) % set x limit
    set(handles.textA, 'String', num2str(negA),'FontSize',10); % display A
    set(handles.textB, 'String', num2str(negB),'FontSize',10); % display B
end

title(foldDir(filePointer).name); % display particle name as title
set(handles.textSpeed, 'String', num2str(partSpeed{filePointer})) % display particle speed
ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % set ylim to nominal baseline and just above max Y
set(handles.textDateTime, 'String', datestr(partTime{filePointer})) % display particle

function resetView_Callback(hObject, eventdata, handles) % reset view of plot normal
global filePointer Data posMZ negMZ polarity NumPoints % declare global variables
if polarity{filePointer} == 0 % determine if particle is positive or negative
	xlim([0 posMZ(NumPoints)]); % reset to posMZ upper limit
else
    xlim([0 negMZ(NumPoints)]); % reset to negMZ upper limit
end
ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % reset y-axis

function zoomToggle_Callback(hObject, eventdata, handles) % toggle zoom for plot
g = get(gca); % get axes handle
outZoom = zoom(gcf); % determine zoom state
if strcmp('off',outZoom.Enable) == 1 % toggle on
    zoom on
else % toggle off
    zoom off
    xlim(g.XLim); % attempt to get around reset of view after exiting zoom--still unsuceesful
    ylim(g.YLim); % attempt to get around reset of view after exiting zoom--still unsuceesful
end

function setX_Callback(hObject, eventdata, handles) % set X axis limts manually
temp_min = inputdlg('Enter minimum X-value'); % min value input
temp_max = inputdlg('Enter maximum X-value'); % max value input 
xlim([str2num(char(temp_min)) str2num(char(temp_max))]) % set limits from user inputs

function setY_Callback(hObject, eventdata, handles) % set Y axis limts manually
temp_min = inputdlg('Enter minimum Y-value'); % min value input
temp_max = inputdlg('Enter maximum Y-value'); % max value input
ylim([str2num(char(temp_min)) str2num(char(temp_max))]) % set limits from user inputs

function setCalibrationValues_Callback(hObject, eventdata, handles) % set A and B calibration values manually
global filePointer Data foldDir posMZ negMZ negA negB posA posB polarity NumPoints % declare global variables
tempA = inputdlg('Enter new A value'); % user input for A
tempB = inputdlg('Enter new B value'); % user input for B
if polarity{filePointer} == 0 % positive
    posA = str2num(char(tempA)); % reset A
    posB = str2num(char(tempB)); % reset B
    posMZ = (posA*([1:NumPoints]-posB).^2)'; % recalculate posMZ array
    set(handles.textA, 'String', num2str(posA),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(posB),'FontSize',10); % display new B
    plot(posMZ,Data{filePointer}); % replot spectrum
    xlim([0 posMZ(NumPoints)]); % set xlim
else % negative
    negA = str2num(char(tempA)); % user input for A
    negB = str2num(char(tempB)); % user input for B
    negMZ = (negA*([1:NumPoints]-negB).^2)'; % recalculate negMZ array
    set(handles.textA, 'String', num2str(negA),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(negB),'FontSize',10); % display new B
    plot(negMZ,Data{filePointer}); % replot spectrum
    xlim([0 negMZ(NumPoints)]); % set xlim
end
title(foldDir(filePointer).name); % display particle name as title
ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % set ylim to nominal baseline and just above max Y

function addPoint_Callback(hObject, eventdata, handles) % add point for calibration--very important function
global filePointer posMZ negMZ polarity particleString particleSaveString % declare global variables
[X,~] = ginput(1); % activate user input through mouse
if polarity{filePointer} == 0
    tmp = abs(posMZ-X);
    [~, idx] = min(tmp); %index of closest value
else
    tmp = abs(negMZ-X);
    [~, idx] = min(tmp); %index of closest value
end
set(handles.textPoints, 'String', sprintf('%s, %s',num2str(X),num2str(idx))); % output index and value to textbox
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

function exitFigure_Callback(hObject, eventdata, handles) % quit program
clear filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints  particleString particleSaveString exportPath % clear global variables
close(gcf) % close figure

function uitable1_CellEditCallback(hObject, eventdata, handles) % make sure particle data is recorded when editing calibration table
global filePointer particleString particleSaveString % declare global variables
[X] = eventdata.Indices; % get new value put int
particleSaveString{X(1)} = particleString{filePointer}; % get string for current particle and save it

function getPoint_Callback(hObject, eventdata, handles) % get point
global polarity filePointer posMZ negMZ
[X,~] = ginput(1); % user input prompt
if polarity{filePointer} == 0 % positive
    tmp = abs(posMZ-X);
    [~, idx] = min(tmp); %index of closest value
else % negative
    tmp = abs(negMZ-X);
    [~, idx] = min(tmp); %index of closest value
end
set(handles.textPoints, 'String', sprintf('%s, %s',num2str(X),num2str(idx))); % display points in textbox

function loadCalibration_Callback(hObject, eventdata, handles)
global filePointer foldDir defaultPath Data posMZ negMZ negA negB posA posB polarity NumPoints % clear global variables
[filename, pathName,~] = uigetfile('*.cal','load calibration data from .cal',defaultPath) % user selects folder to load
fid = fopen([pathName filename ]); % load config values
defaultPath = fgetl(fid); % load Path 
posA = str2num(fgetl(fid)); % load cal value
posB = str2num(fgetl(fid)); % load cal value
negA = str2num(fgetl(fid)); % load cal value
negB = str2num(fgetl(fid)); % load cal value
fclose(fid);
NumPoints = 15000; % declare constant
posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg
if polarity{filePointer} == 0 % positive
    set(handles.textA, 'String', num2str(posA),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(posB),'FontSize',10); % display new B
    plot(posMZ,Data{filePointer}); % replot spectrum
    xlim([0 posMZ(NumPoints)]); % set xlim
else % negative
    negMZ = (negA*([1:NumPoints]-negB).^2)'; % recalculate negMZ array
    set(handles.textA, 'String', num2str(negA),'FontSize',10); % display new A
    set(handles.textB, 'String', num2str(negB),'FontSize',10); % display new B
    plot(negMZ,Data{filePointer}); % replot spectrum
    xlim([0 negMZ(NumPoints)]); % set xlim
end
title(foldDir(filePointer).name); % display particle name as title
ylim([median(Data{filePointer}), max(Data{filePointer})*1.05]); % set ylim to nominal baseline and just above max Y

function saveCalibration_Callback(hObject, eventdata, handles)
global defaultPath negA negB posA posB % clear global variables
[filename, pathName,~] = uiputfile('*.cal','save calibration data as .cal',defaultPath); % user selects folder to load

fid = fopen([pathName filename],'w') %load file for writing
outFormat = '%s \n'; % output format
fprintf(fid,outFormat,defaultPath); % write Path
outFormat = '%d \n'; % output format
fprintf(fid,outFormat,posA); % write cal value
fprintf(fid,outFormat,posB); % write cal value
fprintf(fid,outFormat,negA); % write cal value
fprintf(fid,outFormat,negB); % write cal value
fclose(fid);

function figure1_DeleteFcn(hObject, eventdata, handles)
clear filePointer defaultPath Data partSpeed partTime posMZ negMZ foldDir negA negB posA posB polarity NumPoints...
    particleString particleSaveString exportPath % clear global variables

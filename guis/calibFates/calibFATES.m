function calibFATES
% variable declaration
    NumPoints = 15000; % constant for number of points in spectrum
    calibToggle = 0; % calibration toggle set to default
    
    posA = 1.88E-6; % set default values for calibration parameters
    posB = 0; % set default values for calibration parameters
    negA = 1.4E-6; % set default values for calibration parameters
    negB = 0; % set default values for calibration parameters
    
    posMZ = 1:NumPoints; % uncalibrated spectrum 
    negMZ = 1:NumPoints; % uncalibrated spectrum
    for j = 1:20
        particleSaveString{j} = {''}; % create a cell of empty strings to be overwritten during the recording of which calibration points came from which particle
    end
    particleFilePath = [];
    filePointer = [];
    Data = [];
    partSpeed = [];
    partTime = [];
    foldDir = [];
    polarity = [];
    filename = [];
    spectraIdx = [];
	loadSpectra;% load spectrum from user input%     
   
    % initialize figure
    gf = figure('Visible','off');
    set(gf,'units','normalized');
    set(gf,'position',[0 0 1 1]);
    set(gf,'KeyPressFcn',@keyPressParse);

    topMenu = uimenu(gf,'Label','calibFATES');
    mOpen = uimenu(topMenu,'Label','Open spectrum','Accelerator','u','Callback',@menuCallback);
    mCal = uimenu(topMenu,'Label','Calibration');
    mAddPoint = uimenu(mCal,'Label','Add calibration point','Accelerator','d','Callback',@menuCallback);
    mSaveCal = uimenu(mCal,'Label','Save current calibration values','Accelerator','e','Callback',@menuCallback);
	mLoadCal = uimenu(mCal,'Label','Load calibration values','Accelerator','l','Callback',@menuCallback);
    mSetCal = uimenu(mCal,'Label','Manually set calibration values','Accelerator','b','Callback',@menuCallback);
    mView = uimenu(topMenu,'Label','View');
    mResetView = uimenu(mView,'Label','Reset view','Accelerator','r','Callback',@menuCallback);
    mSetX = uimenu(mView,'Label','Set X axis','Accelerator','t','Callback',@menuCallback);
    mSetY = uimenu(mView,'Label','Set Y axis','Accelerator','y','Callback',@menuCallback);

    posCheckLabel = uicontrol('Parent',gf,'Style','text','String','Positive Spectra','HorizontalAlignment','left','Units','normalized','Position',[0.795,.105,0.065,.02],'FontSize',11);
    posCheck = uicontrol('Parent',gf,'Style','checkbox','Units','normalized','Position',[0.78,0.105,0.01,.02],'Value',1,'Callback',@boxCheck);
    negCheckLabel = uicontrol('Parent',gf,'Style','text','String','Negative Spectra','HorizontalAlignment','left','Units','normalized','Position',[0.795,.08,0.065,.02],'FontSize',11);
    negCheck = uicontrol('Parent',gf,'Style','checkbox','Units','normalized','Position',[0.78,0.08,0.01,.02],'Value',1,'Callback',@boxCheck);
    
    tgCalButton = uicontrol('Parent',gf,'Style','pushbutton','String','Toggle Calibration','Units','normalized','Position',[0.78,0.33,0.15,0.05],'FontWeight','Bold','Callback',@toggleCal,'KeyPressFcn',@keyPressParse);
    clearTableButton = uicontrol('Parent',gf,'Style','pushbutton','String','Clear Table','Units','normalized','Position',[0.78,0.28,0.07,0.04],'Callback',@clearTable,'KeyPressFcn',@keyPressParse);
    calButton = uicontrol('Parent',gf,'Style','pushbutton','String','Calibrate','Units','normalized','Position',[0.86,0.28,0.07,0.04],'Callback',@calibrateButtonPress,'KeyPressFcn',@keyPressParse);
    clearPlotButton = uicontrol('Parent',gf,'Style','pushbutton','String','Clear Plot','Units','normalized','Position',[0.78,0.235,0.07,0.04],'Callback',@clearPlot,'KeyPressFcn',@keyPressParse);
    finalCalButton = uicontrol('Parent',gf,'Style','pushbutton','String','Finalize Calibrate','Units','normalized','Position',[0.86,0.235,0.07,0.04],'KeyPressFcn',@saveCalibrationButton);
    
    aLabel = uicontrol('Parent',gf,'Style','text','String','A','Units','normalized','Position',[0.1,.105,0.065,.05],'FontSize',18,'FontWeight','Bold');
    bLabel = uicontrol('Parent',gf,'Style','text','String','B','Units','normalized','Position',[0.16,.105,0.065,.05],'FontSize',18,'FontWeight','Bold');
    speedLabel = uicontrol('Parent',gf,'Style','text','String','Speed','Units','normalized','Position',[0.22,.105,0.065,.05],'FontSize',18,'FontWeight','Bold');
    dateLabel = uicontrol('Parent',gf,'Style','text','String','Date and Time','Units','normalized','Position',[0.32,.105,0.125,.05],'FontSize',18,'FontWeight','Bold');
    partCoordLabel = uicontrol('Parent',gf,'Style','text','Units','normalized','Position',[0.55,0.08,.18,.08],'BackgroundColor','white','FontSize',16);
    
    if polarity{filePointer} == 0;
        aText = uicontrol('Parent',gf,'Style','text','String',num2str(posA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
        bText = uicontrol('Parent',gf,'Style','text','String',num2str(posB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
    else
        aText = uicontrol('Parent',gf,'Style','text','String',num2str(negA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
        bText = uicontrol('Parent',gf,'Style','text','String',num2str(negB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
    end
    
    speedText = uicontrol('Parent',gf,'Style','text','String',num2str(partSpeed{filePointer}),'Units','normalized','Position',[0.22,.05,0.065,.05],'FontSize',12);
    dateText = uicontrol('Parent',gf,'Style','text','String',datestr(partTime{filePointer}),'Units','normalized','Position',[0.32,.05,0.125,.05],'FontSize',12);
   
    getSpectraIdx; % determine list of spectra to toggle from--depends upon positive and negative checkboxes
    
    ax = axes('Parent',gf,'Units','normalized','Position',[0.10 0.2 0.65 0.7]);
    mainPlot = plot(ax,posMZ,Data{filePointer});
    mainPlot.ButtonDownFcn = @printPosition;
    ax.ButtonDownFcn = @printPosition;
    if polarity{filePointer} == 0 
        ax.XLim = [posMZ([1 NumPoints])];
    else
        ax.XLim = [posMZ([1 NumPoints])];
    end
    h = zoom(gf);
    title(filename);

    calTable = uitable(gf); % Have to call and list properties separtely because trying to do it in one line calls 'a' as a uitable peer
    set(calTable,'Data',[zeros(20,2)],'ColumnName',{'Timepoint','MZ'},'Units','normalized','ColumnEditable',[true,true],'Position',[0.78,0.4,0.15,0.50],'CellEditCallback',@editTable);

    hManager = uigetmodemanager(gf);
    
    gf.Visible = 'on';
function loadSpectra % load spectrum from user input

    [filename, pathName,fileIDX] = uigetfile({'*.ams','Open AMS file';'*.amz','Open AMZ file'},'Open AMS file'); % user selects folder to load
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
        %% Put in code to deal with other file format
    end
    foldDir = dir([pathName '*.ams']); % select ams files to read in
    if isempty(foldDir)
        foldDir = dir([pathName '*.amz']); % if no ams are present, check for amz
    end
    filePointer = 1; % create filePointer for keeping track of spectra 
    
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
        particleFilePath{i} = [pathName foldDir(i).name]; % find full file path of each particle for calibration tracking
    end
end

    function saveCalibrationButton(source,eventdata)
        saveCalibration
    end

    function getSpectraIdx  % determine index of particles for flipping through
        polIdx = cell2mat(polarity); % convert cell to matrix for easy sorting
        spectraIdx = []; % declare blank variable
        if posCheck.Value == 1; % if positive spectra box is checked
            spectraIdx = sort([spectraIdx find(polIdx == 0)]);
        end
        if negCheck.Value == 1; % if negative spectra box is checked
            spectraIdx = sort([spectraIdx find(polIdx == 1)]);
        end
    end

    function clearPlot(source,callbackdata)
        ax.NextPlot = 'replace';
        if polarity{filePointer} == 0
            mainPlot = plot(ax,posMZ,Data{filePointer});
            mainPlot.ButtonDownFcn = @printPosition;
            ax.ButtonDownFcn = @printPosition;
            ax.XLim = [posMZ([1 NumPoints])];
        else
            mainPlot = plot(ax,negMZ,Data{filePointer});
            mainPlot.ButtonDownFcn = @printPosition;
            ax.ButtonDownFcn = @printPosition;
            ax.XLim = [negMZ([1 NumPoints])];
        end
        title(foldDir(filePointer).name);
    end

    function toggleCal(source,callbackdata)
         currentXLim = ax.XLim;
        if polarity{filePointer} == 0
            tmpMin = abs(posMZ-currentXLim(1)); % find closest value for min
            tmpMax = abs(posMZ-currentXLim(2)); % find closest value for max
            [~, newXMin] = min(tmpMin); %index of closest value
            [~, newXMax] = min(tmpMax); %index of closest value
        else 
            tmpMin = abs(negMZ-currentXLim(1)); % find closest value for min
            tmpMax = abs(negMZ-currentXLim(2)); % find closest value for max
            [~, newXMin] = min(tmpMin); %index of closest value
            [~, newXMax] = min(tmpMax); %index of closest value
        end
        if calibToggle == 0
            calculateMZ;
            calibToggle = 1;
        else
            posMZ = 1:NumPoints;
            negMZ = 1:NumPoints;
            calibToggle = 0;
        end
        if polarity{filePointer} == 0;
            mainPlot.XData = posMZ;
            ax.XLim = [posMZ(newXMin) posMZ(newXMax)];
        else
            mainPlot.XData = negMZ;
            ax.XLim = [negMZ(newXMin) negMZ(newXMax)];
        end
    end

    function navigateSpectra(direction)% change plot to display different spectra
        if length(filePointer == 0)
            for i = 1:length(foldDir)
                if strcmp(foldDir(i).name,ax.Title) == 1
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
        currentXLim = ax.XLim;
        if polarity{filePointer} == 0
            tmpMin = abs(posMZ-currentXLim(1)); % find closest value for min
            tmpMax = abs(posMZ-currentXLim(2)); % find closest value for max
            [~, newXMin] = min(tmpMin); %index of closest value
            [~, newXMax] = min(tmpMax); %index of closest value
        else 
            tmpMin = abs(negMZ-currentXLim(1)); % find closest value for min
            tmpMax = abs(negMZ-currentXLim(2)); % find closest value for max
            [~, newXMin] = min(tmpMin); %index of closest value
            [~, newXMax] = min(tmpMax); %index of closest value
        end
        if strcmp(direction,'right')        
            if length(filePointer == 1)
                if filePointer == spectraIdx(end) % if at end of folder, don't redraw plot
                else % else redraw plot
                    filePointer = spectraIdx(currentIdx+1); % go to next spectra
                    if polarity{filePointer} == 0
                        mainPlot.XData = posMZ;
                        ax.XLim = [posMZ(newXMin) posMZ(newXMax)];
                        aText.String = num2str(posA);
                        bText.String = num2str(posB);
                    else 
                        mainPlot.XData = negMZ;
                        ax.XLim = [negMZ(newXMin) negMZ(newXMax)];
                        aText.String = num2str(negA);
                        bText.String = num2str(negB);
                    end
                    mainPlot.YData = Data{filePointer};
                    speedText.String = num2str(partSpeed{filePointer});
                    dateText.String = datestr(partTime{filePointer});
                    title(foldDir(filePointer).name);
                end
            end
        elseif strcmp(direction,'left')
            if filePointer == spectraIdx(1) % if at end of folder, don't redraw plot
            else % else redraw plot
                filePointer = spectraIdx(currentIdx-1); % go to next spectra
                if polarity{filePointer} == 0
                    mainPlot.XData = posMZ;
                    ax.XLim = [posMZ(newXMin) posMZ(newXMax)];
                    aText.String = num2str(posA);
                    bText.String = num2str(posB);
                else 
                    mainPlot.XData = negMZ;
                    ax.XLim = [negMZ(newXMin) negMZ(newXMax)];
                    aText.String = num2str(negA);
                    bText.String = num2str(negB);
                end
                mainPlot.YData = Data{filePointer};
                speedText.String = num2str(partSpeed{filePointer});
                dateText.String = datestr(partTime{filePointer});
                title(foldDir(filePointer).name);
            end    
        end
        ax.YLimMode = 'auto';
    end

    function menuCallback(source,callbackdata)
        if strcmp(source.Label,'Open spectrum')
            loadSpectra
            mainPlot = plot(ax,posMZ,Data{filePointer});
            mainPlot.ButtonDownFcn = @printPosition;
            ax.ButtonDownFcn = @printPosition;
            title(filename)
            if polarity{filePointer} == 0;
                aText = uicontrol('Parent',gf,'Style','text','String',num2str(posA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
                bText = uicontrol('Parent',gf,'Style','text','String',num2str(posB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
            else
                aText = uicontrol('Parent',gf,'Style','text','String',num2str(negA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
                bText = uicontrol('Parent',gf,'Style','text','String',num2str(negB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
            end
        elseif strcmp(source.Label,'Add calibration point')
            [X,~] = ginput(1);
            addCalPoint(X)
        elseif strcmp(source.Label,'Save current calibration values')
            saveCalibration
        elseif strcmp(source.Label,'Load calibration values')
            loadCalibration
        elseif strcmp(source.Label,'Manually set calibration values')
            setCalibration
        elseif strcmp(source.Label,'Reset view')
            clearPlot(source,callbackdata)
        elseif strcmp(source.Label,'Set X axis')
            setAxis('X');
        elseif strcmp(source.Label,'Set Y axis')
            setAxis('Y');
        end
        
    end
        
    function setAxis(axis)
        if strcmp(axis,'X')
            XMin = str2num(char(inputdlg('Please enter new X min')));
            XMax = str2num(char(inputdlg('Please enter new X max')));
            if ~isempty(XMin)
                if ~isempty(XMax)
                    ax.XLim = [XMin XMax];
                end
            end
        elseif strcmp(axis,'Y')
            YMin = str2num(char(inputdlg('Please enter new Y min')));
            YMax = str2num(char(inputdlg('Please enter new Y max')));
            if ~isempty(YMin)
                1
                if ~isempty(YMax)
                    ax.YLim = [YMin YMax];
                end
            end
        end
    end

    function calculateMZ % calculate MZ variables for x-axis
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
        posMZ = (posA*([1:NumPoints]-posB).^2)'; % create M/Z for pos
        negMZ = (negA*([1:NumPoints]-negB).^2)'; % create M/Z for neg 
    end

    function clearTable(source,callbackdata)
        calTable.Data = zeros(20,2);
    end

    function editTable(hObject, eventdata) 
        [X] = eventdata.Indices; % get new value put int
        particleSaveString{X(1)} = particleFilePath{filePointer}; % get string for current particle and save it
    end

    function keyPressParse(source,eventData)
        if strcmp(eventData.Key,'m')
            if strcmp(eventData.Modifier,'control')
                if strcmp(h.Enable,'on')
                    h.Enable = 'off';
                elseif strcmp(h.Enable,'off')
                    h.Enable = 'on'
                    [hManager.WindowListenerHandles.Enabled] = deal(false) % HG2
                    gf.KeyPressFcn = @keyPressParse;
                end
            end
        end
        if strcmp(eventData.Key,'leftarrow')
            navigateSpectra('left');
        elseif strcmp(eventData.Key,'rightarrow');
            navigateSpectra('right');
        end
    end
        
    function printPosition(source,callbackdata)
        if callbackdata.Button == 1
            partCoordLabel.String = num2str(callbackdata.IntersectionPoint(1:2));
        elseif callbackdata.Button == 3
            partCoordLabel.String = num2str(callbackdata.IntersectionPoint(1:2));
            X = callbackdata.IntersectionPoint(1);
            addCalPoint(X);
        end
    end 

    function calibrateButtonPress(source,callbackdata)
        if calibToggle == 0
        toggleCal;
        end
            calData = calTable.Data; % get data currently in table
            [A,B] = massCal(calData(calData(:,1)>0,1),calData(calData(:,1)>0,2)); % determine A and B from data in table
            tempMZ = A*([1:NumPoints]-B).^2'; % calculate temp MZ
            ax.NextPlot = 'Add' 
                plot(tempMZ,Data{filePointer}) % plot spectra with new MZ
                title(foldDir(filePointer).name); % display particle name as title
            ax.NextPlot = 'Replace'
                aText.String = num2str(A);
                bText.String = num2str(B);
    end

    function addCalPoint(X) % add point for calibration--very important function
        if polarity{filePointer} == 0
            tmp = abs(posMZ-X);
            [~, idx] = min(tmp); %index of closest value
        else
            tmp = abs(negMZ-X);
            [~, idx] = min(tmp); %index of closest value
        end
        Y = inputdlg('Please enter M/Z'); % get user input MZ for point
        if isempty(Y)
            return
        else
            calCounter = find(calTable.Data(:,1) == 0); % find first value with zero for the timepoint value
            calTable.Data(calCounter(1),1) = idx;
            calTable.Data(calCounter(1),2) = str2num(char(Y));
            particleSaveString{calCounter(1)} = particleFilePath{filePointer};
        end
    end

    function saveCalibration
        aValue = str2num(aText.String);
        bValue = str2num(bText.String);
        if polarity{filePointer} == 0
            posA = aValue;
            posB = bValue;
        else
            negA = aValue;
            negB = bValue;
        end
        if calibToggle == 1
            calculateMZ
        end
            
        
        [fName,pName,~] = uiputfile('*.cal','save calibration data'); % user selects folder to load
        fid = fopen([pName fName],'w');
        
        outFormat = '%d \n'; % output format
        fprintf(fid,outFormat,aValue); % write cal value
        fprintf(fid,outFormat,bValue); % write cal value
        outFormat2 = '%d, %d, ';
        outFormat3 = '%s\n';
        calIdx = calTable.Data(:,1)>0;
        if isempty(calTable.Data(calTable.Data(:,1)>0,1))
        else
            for i = 1:length(calIdx)
                fprintf(fid,outFormat2,calTable.Data(calIdx(i),:));
                fprintf(fid,outFormat3,particleSaveString{calIdx(i)});
            end
        end
        fclose(fid);
    end

    function loadCalibration
        [fName,pName,~] = uigetfile('*.cal','load calibration data'); % user selects folder to load
        fid = fopen([pName fName],'r');
        if polarity{filePointer} == 0
            posA = str2num(fgetl(fid));
            posB = str2num(fgetl(fid));
            aText = uicontrol('Parent',gf,'Style','text','String',num2str(posA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
            bText = uicontrol('Parent',gf,'Style','text','String',num2str(posB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
        else
            negA = str2num(fgetl(fid));
            negB = str2num(fgetl(fid));
            aText = uicontrol('Parent',gf,'Style','text','String',num2str(negA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
            bText = uicontrol('Parent',gf,'Style','text','String',num2str(negB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
        end
        fclose(fid);
    end

    function setCalibration
        tempA = inputdlg('Enter new A value'); % user input for A
        tempB = inputdlg('Enter new B value'); % user input for B
        if polarity{filePointer} == 0
            posA = str2num(char(tempA));
            posB = str2num(char(tempB));
            aText = uicontrol('Parent',gf,'Style','text','String',num2str(posA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
            bText = uicontrol('Parent',gf,'Style','text','String',num2str(posB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
        else
            negA = str2num(char(tempA));
            negB = str2num(char(tempB));
            aText = uicontrol('Parent',gf,'Style','text','String',num2str(negA),'Units','normalized','Position',[0.1,.05,0.065,.05],'FontSize',12);
            bText = uicontrol('Parent',gf,'Style','text','String',num2str(negB),'Units','normalized','Position',[0.16,.05,0.065,.05],'FontSize',12);
        end
    end

    function boxCheck(source, eventData)
        getSpectraIdx;
    end

    function [A,B] = massCal(X,Y,plotFlag) 
    % perform mass calibration with data
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
    end
end

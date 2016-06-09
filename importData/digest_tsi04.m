function digest_tsi04(PKLDir,PK2Dir)
% DIGEST_TSI04 creates PK2 files from TSI data files
% Call as DIGEST_TSI04(PKLDir,PK2Dir)
% where PKLDir is the directory of .pkl, .sem, and .inst data files
%       PK2Dir is the directory to write .pk2 files
% All files in PKLDir and its subdirectories will be processed.
%
% Data directories can contain 1 instrument file for the entire directory, 
% or 1 instrument file for each .pkl file.  If there is one instrument 
% file in a directory, the data are copied to every .pk2 and the base file 
% name can be anything.  If there one instrument file for each .pkl, the 
% base file name must match the .pkl file.
% 
% This function digests raw data files created by the 2004 version of TSI 
% data acquision and analysis software using Matlab functions.
% 
% AMZData file names in SET and PKL files must be sorted by time and 
% have the same order.  
%
% See also READ_INST, READ_SEM, READ_SET, READ_PKL, WRITE_PK2, DIGEST_PK2

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Based on digest_tsi00
% Replaced Perl function (tsi00.pl) with Matlab functions to read and 
% write text data files
% Jonathan O. Allen  2005-01-05

% Converted timer count data into times (s)
% Added particle velocity

% Added warning if more than one .inst file in data directory
% JOA  2007-02-27

% Replace calls to isodatestr with sortrow of datevec for better performance
% Revised for isodatestr return of char matrix instead of cell string vector
% Unwraped for loop over peaks for better performance
% JOA  2008-01-24

% Major revision to read large PKL files
%   Read all SET data
%   Loop over SET data to create PK2 with ~MaxSpec hit particles
%     Read matching missed particle data from SEM
%     Read matching hit particle data from PKL
%     Write PK2
% JOA  2008-03-05

% Replaced textread with textscan for better performance
% JOA  2008-05-21

% Omit leading '+' and insignificant zeros in PK2; reduces size by ~10%
% JOA  2008-07-22

% JC 2013 legacy update

FiveMin = 5/24/60;
OneSec = 1/60/60/24;

% hardware parameters
% these should match values in *.inst
ScatterLength = 0.06; % distance between timing lasers (m)
TimerResolution = 50e-9; % time for one count (s)

if ~isword(PKLDir)
  error('Expecting word for PKLDir');
end
if ~exist(PKLDir,'dir')
  error('Unknow PKLDir %s',PKLDir);
end
if ~isword(PK2Dir)
  error('Expecting word for PK2Dir');
end
if ~exist(PK2Dir,'dir')
  error('Unknow PK2Dir %s',PK2Dir);
end

% maximum number of spectra written to each PK2 file
MaxSpec = 1000;

% datadef text is copied to start of every PK2 file
% revise this to match data structure
i = 0;
i = i + 1;
% list of ytables
DataDefText{i} = 'Table: INST PART SPEC PEAK';
i = i + 1;
% list of fields in INST
DataDefText{i} = 'INST: AvgLaserPower BusyTimeFunction BusyTimeParam DaCalibFunction DaCalibParam InstCode InstDesc InstName ExpDesc ExpName PreProcDesc MinHeight MinArea OpName OpDesc SampleFlow ScatterLength TimerResolution';
i = i + 1;
% list of fields in PART
DataDefText{i} = 'PART: Time Velocity ScatterDelay BusyTime PositionInFolder FastScatter LaserPower';
i = i + 1;
% list of fields in SPEC
DataDefText{i} = 'SPEC: Polarity';
i = i + 1;
% list of fields in PEAK
DataDefText{i} = 'PEAK: MZ Area Height BlowScale';

% construct list of related text data files
% DataFile is structure with fields for related text data files 
% The fields are: Inst, SEM, SEF, SET, PKL 
DataFile = struct([]);
DataFile = update_datafile(PKLDir,DataFile);

for i = 1:length(DataFile)
  %disp(DataFile(i).Inst);
  %disp(DataFile(i).SEM);
  %disp(DataFile(i).SEF);
  %disp(DataFile(i).SET);
  %disp(DataFile(i).PKL);
  %disp(' ');

  [SETData,AMZNameSET] = read_set(DataFile(i).SET);
  SETData = [SETData(:,2) SETData(:,3) SETData(:,4) SETData(:,1) SETData(:,5) SETData(:,6)];

  % round second to avoid issues related to sorting floats
  %SETTimeVec = datevec(SETData(:,1));
  %SETTimeVec(:,end) = round(SETTimeVec(:,end));
  %[SETTimeVec,idx] = sortrows(SETTimeVec);
  %AMZNameSET = AMZNameSET(idx);

  % loop over groups of MaxSpec hit particles using SET data
  NumSpec = size(SETData,1);
  NumGroups = ceil(NumSpec/MaxSpec);
  SETStartIdx = zeros(NumGroups,1);
  SETStopIdx  = zeros(NumGroups,1);
  SETStartTime = zeros(NumGroups,1);
  SETStopTime  = zeros(NumGroups,1);
  AMZNameSET2 = cell(NumGroups,1);

  j = 0; % group index
  si1 = 0; % spectra index
  while si1 < NumSpec
    j = j + 1;
    si0 = si1 + 1;
    SETStartIdx(j) = si0;
    if j == 1;
      SETStartTime(j) = SETData(1,1) - FiveMin;
    else
      SETStartTime(j) = SETStopTime(j-1) + OneSec;
    end      
    si1 = si1 + MaxSpec;
    if si1+1 > NumSpec
      si1 = NumSpec;
    else
      % include all spectra in same second in group
      while SETData(si1,1) > SETData(si1+1,1) - OneSec/2 & SETData(si1,1) < SETData(si1+1,1) + OneSec/2 
        si1 = si1 + 1;
        if si1+1 > NumSpec
          si1 = NumSpec;
          break;
        end
      end
    end
    SETStopIdx(j) = si1;
    SETStopTime(j) = SETData(si1,1);
  end

  for j = 1:length(SETStartIdx)
    if SETStartIdx(j) == 0
      % empty last group
      continue;
    end
    
    NumSET = SETStopIdx(j) - SETStartIdx(j) + 1;
    
    % get missed particle data
    % SE*Data matrices have rows with format:
    % [Count Time ScatterDelay BusyTime FastScatter LaserPower]
    SEMData = read_sem(DataFile(i).SEM,SETStartTime(j),SETStopTime(j),0);
    if ~isempty(DataFile(i).SEF)
      SEFData = read_sem(DataFile(i).SEF,SETStartTime(j),SETStopTime(j),1);
    else
      SEFData = [];
    end
    
    % move Count column so column order is
    % [Time ScatterDelay BusyTime Count FastScatter LaserPower]
    % sort data on Time  
    SEMData = [SEMData(:,2) SEMData(:,3) SEMData(:,4) SEMData(:,1) SEMData(:,5) SEMData(:,6)];
    % round second to avoid issues related to sorting floats
    %SEMTimeVec = datevec(SEMData(:,1));
    %SEMTimeVec(:,end) = round(SEMTimeVec(:,end));
    %[SEMTimeVec,idx] = sortrows(SEMTimeVec);
    %SEMData = SEMData(idx,:);

    if ~isempty(SEFData)
      SEFData = [SEFData(:,2) SEFData(:,3) SEFData(:,4) SEFData(:,1) SEFData(:,5) SEFData(:,6)];
      % round second to avoid issues related to sorting floats
      %SEFTimeVec = datevec(SEFData(:,1));
      %SEFTimeVec(:,end) = round(SEFTimeVec(:,end));
      %[SEFTimeVec,idx] = sortrows(SEFTimeVec);
      %SEFData = SEFData(idx,:);
    end

    % PKLData rows correspond to peaks; columns are for the data fields
    % in this format:
    %   [ParticleSerial Polarity MassToCharge PeakArea PeakHeight BlowScale]
    [PKLData,AMZNamePKL] = read_pkl(DataFile(i).PKL,AMZNameSET(SETStartIdx(j):SETStopIdx(j)));

    % match peak data
    SETPeak1 = zeros(NumSET,1);
    SETPeak2 = zeros(NumSET,1);

    [UPKLSerial,UPKLSerialIdx2] = unique(PKLData(:,1),'rows','legacy'); % JC 2013
    UPKLSerialIdx1 = [1; UPKLSerialIdx2(1:end-1)+1];

    l = 1;
    for k = 1:NumSET
      if UPKLSerial(l) == (k + SETStartIdx(j) - 1)
        SETPeak1(k) = UPKLSerialIdx1(l);
        SETPeak2(k) = UPKLSerialIdx2(l);
        l = l + 1;
      end
    end
    
    % get instrument text for transcription and InstCode for PK2 file name
    InstCode = '';

    fid = fopen(DataFile(i).Inst,'r');
    if fid == -1
      error('Cannot find Inst file %s',DataFile(i).Inst);
    end

    InstText = textscan(fid,'%s','delimiter','\n');
    fclose(fid);

    % remove blank lines
    InstText = InstText{:};
    InstText_empty = cellfun('isempty',InstText);
    InstText = InstText(~InstText_empty);
    
    for k = 1:length(InstText)
      tline = InstText{k};
      if ~isempty(tline)
        [firstword,tline] = strtok(tline);
        if ~isempty(regexp(upper(firstword),'^INSTCODE', 'once' ));
          % remove '='
          [junk,InstCode] = strtok(tline);
          InstCode = trim(InstCode);
        end
      end
    end
    if isempty(InstCode)
      error('InstCode not found in %s',DataFile(i).Inst);
    end  

    % create file name
    dstr = datestr(SETData(SETStartIdx(j),1),30);
    PK2File = fullfile(PK2Dir,sprintf('%s%s.pk2',InstCode,dstr));
    PK2fid = fopen(PK2File,'w');

    disp(sprintf('Creating %s',PK2File));
    disp([datestr(SETStartTime(j)) ' - ' datestr(SETStopTime(j))]);
    
    % write creation comments in form:
    % Created by DIGEST_TSI04.m on 03-Jan-2003 19:42:39 from source files:
    %   c:\data\biis\raw\Jan12\Bakersfield.inst
    %   c:\data\biis\raw\Jan12\zzzb1.pkl
    %   c:\data\biis\raw\Jan12\zzzb1.sem
    %   [no .sef]
    fprintf(PK2fid,'%% Created by DIGEST_TSI04.m on %s from source files:\n',isodatestr(now));
    fprintf(PK2fid,'%%   %s\n', DataFile(i).Inst);
    fprintf(PK2fid,'%%   %s\n', DataFile(i).SEM);
    if ~isempty(DataFile(i).SEF)
      fprintf(PK2fid,'%%   %s\n', DataFile(i).SEF);
    end
    fprintf(PK2fid,'%%   %s\n', DataFile(i).SET);
    fprintf(PK2fid,'%%   %s\n\n', DataFile(i).PKL);
    
    % write data definition
    fprintf(PK2fid,'%% Data Definition\n');
    fprintf(PK2fid,'%s\n', DataDefText{:});
    fprintf(PK2fid,'\n');
    
    % transcribe INST data
    fprintf(PK2fid,'%s\n', InstText{:});
    fprintf(PK2fid,'\n');

    % collect particle data in time order
    if ~isempty(SEFData)
      AllData = [SETData(SETStartIdx(j):SETStopIdx(j),:); SEMData; SEFData];
    else
      AllData = [SETData(SETStartIdx(j):SETStopIdx(j),:); SEMData];
    end
    [AllData,AllIdx] = sortrows(AllData);

    % track rows of SETData in AllData
    SETIdx2 = zeros(size(AllData,1),1);
    SETIdx2(1:NumSET) = 1:NumSET;
    SETIdx2 = SETIdx2(AllIdx);
    
    TimeText = isodatestr(AllData(:,1));

    % convert timer counts (50 ns) to time (s) for ScatterDelay and BusyTime
    ScatterDelay = AllData(:,2) * TimerResolution;
    BusyTime = AllData(:,3) * TimerResolution;

    % create velocity data
    Velocity = ScatterLength ./ ScatterDelay;
    
    for k = 1:size(AllData,1)

      % write particle data in format:
      % Time Velocity ScatterDelay BusyTime PositionInFolder FastScatter LaserPower
      fprintf(PK2fid,'%s %g %g %g %i %i %g\n',TimeText(k,:),Velocity(k),ScatterDelay(k),BusyTime(k),AllData(k,4:6));
      
      if SETIdx2(k)
        % hit particle
        % find matching PKLData
        % PKLData in format:
        % [ParticleSerial Polarity MassToCharge PeakArea PeakHeight BlowScale]

        %        PKLPeakIdx = binary_search(PKLData(:,1),'==',SETIdx2(k));
        % disp([SETPeak1(SETIdx2(k)) SETPeak2(SETIdx2(k))]);
        PKLPeakIdx = SETPeak1(SETIdx2(k)):SETPeak2(SETIdx2(k));
        NegPeakIdx = PKLPeakIdx(find(PKLData(PKLPeakIdx,2) == 0));
        PosPeakIdx = setxor(PKLPeakIdx,NegPeakIdx,'rows','legacy');% JC 2013
        
        if ~isempty(PosPeakIdx)
          % write spectrum data in format:
          % > Polarity
          fprintf(PK2fid,'> 1\n');
          
          % write peak data in format:
          % >> MZ Area Height BlowScale
          fprintf(PK2fid,'>> %g %i %i %i\n',PKLData(PosPeakIdx,3:6)');
        end

        if ~isempty(NegPeakIdx)
          % write spectrum data in format:
          % > Polarity
          fprintf(PK2fid,'> 0\n');
          
          % write peak data in format:
          % >> MZ Area Height BlowScale
          fprintf(PK2fid,'>> %g %i %i %i\n',PKLData(NegPeakIdx,3:6)');
        end
      end
    end

    fclose(PK2fid);
  end
end

return
    


function DataFile = update_datafile(PKLDir,DataFile)
% construct list of related text data files
% DataFile is structure with fields for related text data files 
% The fields are: Inst, SEM, SEF, SET, PKL 

DirList = sort_struct(dir(PKLDir));

% find text data files in PKLDir and subdirectories
InstFileCount = 0;
InstFile = '';
for i = 1:length(DirList)

  if DirList(i).isdir
    
    if strmatch(DirList(i).name,'.')
      % find instrument files in current directory
      for j = 1:length(DirList)
        % check for .inst files
        FileName = fullfile(PKLDir,DirList(j).name);
        [PathName BaseFile Extension] = fileparts(FileName);
        
        if strcmp(upper(Extension),'.INST')
          InstFileCount = InstFileCount + 1;
          InstFile = FileName;
          if InstFileCount > 1
            warning('More than one .inst file in data directory %s',PKLDir);
            break;
          end
        end
      end

    elseif strmatch(DirList(i).name,'..')
      % ignore 
    else
      % process files in subdirectories
      DataFile = update_datafile(fullfile(PKLDir,DirList(i).name),DataFile);
    end
    
  else
    % check for .pkl data files
    PKLFile = fullfile(PKLDir,DirList(i).name);
    [PathName BaseFile Extension] = fileparts(PKLFile);

    if strcmp(upper(Extension),'.PKL')
      Missing = 0;

      % find related files
      if InstFileCount == 0
        disp(sprintf('Cannot digest %s; missing .inst file',PKLFile));
        Missing = 1;
        
      elseif InstFileCount > 1
        % look for new .inst only if more than one in directory
        InstFile = fullfile(PathName,[BaseFile '.inst']); 
        if ~exist(InstFile,'file') 
          disp(sprintf('Cannot digest %s; missing .inst file',PKLFile));
	  Missing = 1;
        end
      end
      
      SEMFile = fullfile(PathName,[BaseFile '.sem']); 
      if ~exist(SEMFile,'file') 
        disp(sprintf('Cannot digest %s; missing .sem file',PKLFile));
	Missing = 1;
      end

      SEFFile = fullfile(PathName,[BaseFile '.sef']); 
      if ~exist(SEFFile,'file') 
        % sef is optional
        SEFFile = '';
      end

      SETFile = fullfile(PathName,[BaseFile '.set']); 
      if ~exist(SETFile,'file') 
        disp(sprintf('Cannot digest %s; missing .set file',PKLFile));
	Missing = 1;
      end

      if ~Missing
	% record data file names
	idx = length(DataFile) + 1;
	DataFile(idx).Inst = InstFile;
	DataFile(idx).SEM  = SEMFile;
	DataFile(idx).SEF  = SEFFile;
	DataFile(idx).SET  = SETFile;
	DataFile(idx).PKL  = PKLFile;
      end
    end
  end
end

return

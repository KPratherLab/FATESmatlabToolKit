%read inst file
function [InstID,NewInst] = parse_inst(InstFile)
% PARSE_INST reads text file of instrument data and updates INST
% Call as [InstID,NewInst] = PARSE_INST(DataDef,InstFile)
% InstID    instrument ID for InstFile
% NewInst   true if this instrument is different from the prior one
% InstFile  text file containing instrument data
%
% InstFile contains one row of instrument data.  It has 2 types of lines 
% 1. comments that start in '%'
% 2. data lines with the form 'Field = Value'
%
% Every data folder has an instrument file; since many data folders have the same 
% instrument data, a new instrument is created only when data in the 
% InstFile differs from a previous InstFile.

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  24 Mar 00

% changed InstID setting to use OpCode
% increment InstID for new InstCode-OpCode-StudyCode combinations
% JOA  17 Feb 01

% fixed increment of InstID for new InstCode-OpName-ExpName combinations
% JOA  24 Dec 01

% JC 2013 legacy update
% PFR/CMS 2016 updated for streamlined YAADA

% global STUDY INST
%PFR
%fprintf('INFO, parse_inst, starting in parseinst for %s\n',InstFile);

%add new row to inst
if ~isempty(INST(1).InstID)
    RowIdx = length(INST)+1;  %PFR numrow(INST) + 1;
else
    RowIdx = 1;
end;

%PFR instid is now integeer InstID = instid(STUDY.LastInstID) + 1;
InstID = STUDY.LastInstID + 1;
%fprintf('INFO, parse_inst, InstID set to %i for rowindex: %i\n',InstID,RowIdx);

%PFR replace the col by col field creation with addig a struct element
%INST(RowIdx,'InstID') = InstID;
INST(RowIdx).InstID=InstID;

%read in inst file
fldname_set=fields(INST); %PFR get field names
instText = textscan(fid,'%s','Delimiter','\n'); %read in inst file
instText = instText{1}(~cellfun('isempty',instText{1})); %get rid of empty lines
instText = strtrim(instText); %get rid of trailing and leading white spaces
instText = instText(~strncmp(instText,'%',1)); %get rid of lines that are comments
instText = regexp(instText,'\=','split'); %split string
instText = strtrim(instText); %get rid of any trailing/leading white spaces

%check INST files for field names to be read in and add data to INST
%structure
for i = 1:size(instText,1)
    hasfld=strcmpi(instText{i,1}{1,1},fldname_set); %compare to INST fields
    if any(hasfld)
        useFld = fldname_set(hasfld); %get field name
        %add data to INST structure
        num = str2num(instText{i,1}{1,2}); 
        if isempty(num)
            str2eval=sprintf('INST(RowIdx).%s=''%s'';',useFld{1},instText{i,1}{1,2});
            eval(str2eval)
        else
            str2eval=sprintf('INST(RowIdx).%s=num;',useFld{1});
            eval(str2eval)
        end
    end;
end
fclose(fid);

% reset partid
%PFR INST(RowIdx,'LastPartID') = partid(InstID,0);
%         id2set                 =partid(InstID,0);
INST(RowIdx).LastPartID=0;  %PFR, lastid gets set later
INST(RowIdx).DAcalibFUNCTION = {INST(RowIdx).DAcalibFUNCTION};
INST(RowIdx).DAcalibPARAM = {INST(RowIdx).DAcalibPARAM};
%  disp('INFO, parse_inst, rowidx > 1, so about to check new instrment diff than previous')

% check if new instrument is different from previous instruments
NewInst = 1; %assume new instrument
if RowIdx > 1
%PFR change search table to findmatch structure
%  MatchInst = search(INST.InstCode,'==',char(INST(RowIdx).InstCode));
%  MatchOp   = search(INST.OpName,'==',char(INST(RowIdx).OpName));
%  MatchExp  = search(INST.ExpName,'==',char(INST(RowIdx).ExpName));

%PFR, note we don't put [] around INST.InstCode b/c we want a char array
  MatchInst = strcmpi(INST(RowIdx).InstCODE, {INST(:).InstCODE}); 
  MatchExp  = strcmpi(INST(RowIdx).ExpNAME, {INST(:).ExpNAME});

  Match = MatchInst & MatchExp; % JC 2013
  indexMatch = find(Match);
  if length(indexMatch) > 1
    % inst matches with itself, so test for more than 1 matches
    %fprintf('INFO, parse_inst, match len>1 so newinst set to 0')
    NewInst = 0;
    InstID = INST(indexMatch(1)).InstID; %set InstID
    % INST files can contain different DA parameters for the same experiment
    %(size cals change) Don't want to create new instID for each size cal,
    % so store multiple size cal under single InstID to previous so that 
    % new particles will be calibrated appropriately
    numCal = length(INST(indexMatch(1)).DAcalibFUNCTION);
    %check to see if size cal is different
    matchFunction = ~isequal(INST(RowIdx).DAcalibFUNCTION{1},INST(indexMatch(1)).DAcalibFUNCTION{numCal});
    matchParam = ~isequal(INST(RowIdx).DAcalibPARAM{1},INST(indexMatch(1)).DAcalibPARAM{numCal});
    if matchFunction || matchParam
        INST(indexMatch(1)).DAcalibFUNCTION{numCal+1} = INST(RowIdx).DAcalibFUNCTION;
        INST(indexMatch(1)).DAcalibPARAM{numCal+1} = INST(RowIdx).DAcalibPARAM;
    end
  end
end %end rowidx >1

if NewInst
    % update LastInstID
    %  disp('INFO, parse_inst, new inst')
    STUDY.LastInstID = InstID;
else
    % remove duplicate instrument data
    %  fprintf('INFO, parse_inst, removing duplicate inst data 1:%i -1\n',RowIdx)
    INST = INST(1:RowIdx-1);
    
end

return
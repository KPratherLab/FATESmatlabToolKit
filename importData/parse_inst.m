function [InstID,NewInst] = parse_inst(InstFile)
% PARSE_INST reads file containing instrument and experiment data and updates INST structure
% Call as [InstID,NewInst] = PARSE_INST(InstFile)
% InstID    instrument ID for InstFile
% NewInst   true if this instrument is different from a prior one in the study
% InstFile  text file containing instrument data
%
% Each line of InstFile contains one row of instrument data.  It has 2 types of lines 
% 1. comments that start in '%'
% 2. data lines with the form 'Field = Value'
%
% Every data folder has an instrument file; since many data folders have the same 
% instrument data, a new instrument is created only when an InstFile has a unique
% InstCode and ExpName pair when compared to a previous InstFile.

global STUDY INST

%add new row to inst
if ~isempty(INST(1).InstID)
    RowIdx = length(INST)+1;  
else
    RowIdx = 1;
end;
InstID = STUDY.LastInstID + 1;
INST(RowIdx).InstID=InstID;

%read in inst file
instText = read_inst(InstFile{1});

%check INST files for field names to be read in and add data to INST structure
fldname_set=fieldnames(INST); %get field names
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

% refine data to INST structure
INST(RowIdx).LastPartID=0;  %lastid gets set later
INST(RowIdx).DAcalibFUNCTION = {INST(RowIdx).DAcalibFUNCTION};
INST(RowIdx).DAcalibPARAM = {INST(RowIdx).DAcalibPARAM};

% check if new instrument is different from previous instruments
NewInst = 1; %assume new instrument
if RowIdx > 1
    %note we don't put [] around INST.InstCode b/c we want a char array
    %check to see if inst matches
    MatchInst = strcmpi(INST(RowIdx).InstCODE, {INST(:).InstCODE}); %check if InstCODE matches
    MatchExp  = strcmpi(INST(RowIdx).ExpNAME, {INST(:).ExpNAME}); %check if ExpNAME matches    
    Match = MatchInst & MatchExp; % a match requires both to match
    indexMatch = find(Match); %find row where match occurs
    
    if length(indexMatch) > 1
        % inst matches with itself, so test for more than 1 matches
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
        %add new calibration parameters
        if matchFunction || matchParam
            INST(indexMatch(1)).DAcalibFUNCTION{numCal+1} = INST(RowIdx).DAcalibFUNCTION;
            INST(indexMatch(1)).DAcalibPARAM{numCal+1} = INST(RowIdx).DAcalibPARAM;
        end
    end
    
end %end rowidx >1

if NewInst
    % update LastInstID
    STUDY.LastInstID = InstID;
else
    % remove duplicate instrument data
    INST = INST(1:RowIdx-1);    
end

return
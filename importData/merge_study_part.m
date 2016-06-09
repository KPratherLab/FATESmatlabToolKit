function merge_study_part(inst1_set,inst2_set)
% Call as merge_study_part(inst1_set,inst2_set)
% MERGE_YDB 
%    inst1_set and inst2_set are the study1 inst ids and study2 inst ids that match up,
%   This program will go thru each pair and change inst ids and part ids in study2
%         to follow the last part id for study1 
%   Then it will rewrite all data into merge filename (ie the study3)
%  All instids in study2 need to be assigned to an instid in study1,
%  however study1 can have instids that are not used in reassigning study2
%
global YAADA STUDY DATADEF  runbatch  STUDY1 STUDY2 STUDY3 PARTdataMat PARTidMat PARTdataMat1 PARTidMat1 PARTdataMat2 PARTidMat2 PARTidMiss2 INST INST1 INST2 INST3 PEAK1 PEAK2 PEAK3 PEAK PARTidFlds PARTdataFlds1 PARTdataFlds; 
if nargin > 2
  error('too many arguments, Call as init_study(StudyName1,StudyName2)');
end

if (length(inst1_set) ~= length(inst2_set))
    error('inst1 set and inst2 set must be same length');
end;
%to renum partids we need to match inst ids and get the last part id number
%from study 1

%have to update PEAK info structure, but see below also wehre ids are
%updated
for i=1:length(PEAK2)
    PEAK2(i).rowcnt=PEAK2(i).rowcnt+PEAK1(end).rowcnt;
end;

%get inst and partid data from PEAK
instPEAK1 = [PEAK1.INSTID]';
partidPEAK1 = [PEAK1.PARTID]';
instPEAK2 = [PEAK2.INSTID]';
partidPEAK2 = [PEAK2.PARTID]';
    
%% redo PART data
New_PartIDs = zeros(length(inst1_set),3);
inst2_idx = cell(1,length(inst1_set));
inst2_idxMiss = cell(1,length(inst1_set));

fprintf('INFO,merge study, changing PartIds for study2 ...\n'); 

for i=1:length(inst1_set)  %loop thru inst2 ids for changing partids 
        partid_2start     = partidPEAK1(find(instPEAK1==inst1_set(i),1,'last')); %get last partid for instid in study1 
        inst2_idx{i}        = PARTidMat2(:,PARTidFlds.INSTID)==inst2_set(i); %get all rows from study2 part matrix
        inst2_idxMiss{i}        = PARTidMiss2(:,PARTidFlds.INSTID)==inst2_set(i); %get all rows from study2 part matrix       
        New_PartIDs(i,:) = [inst2_set(i) inst1_set(i) partid_2start]; %save info for PEAK        
        
        %update row conts in the PEAK row,inst/part id structure,
        % this is used sometimes to help find sets of peak rows
        changePEAK2 = find(instPEAK2 == inst2_set(i));
        for ir=1:length(changePEAK2)
            PEAK2(changePEAK2(ir)).INSTID = inst1_set(i);
            PEAK2(changePEAK2(ir)).PARTID = PEAK2(changePEAK2(ir)).PARTID+partid_2start;
        end;       
end;   

%change instid and partid part in study 2
for i = 1:length(inst1_set)
    PARTidMat2(inst2_idx{i},PARTidFlds.PARTID)= PARTidMat2(inst2_idx{i},PARTidFlds.PARTID)+New_PartIDs(i,3); %renumber study2 partids
    PARTidMiss2(inst2_idxMiss{i},PARTidFlds.PARTID)= PARTidMiss2(inst2_idxMiss{i},PARTidFlds.PARTID)+New_PartIDs(i,3); %renumber study2 partids
    PARTidMat2(inst2_idx{i},PARTidFlds.INSTID) = New_PartIDs(i,2); %renumber study2 partids
    PARTidMiss2(inst2_idxMiss{i},PARTidFlds.INSTID)= New_PartIDs(i,2); %renumber study2 partids
end

PARTidMat=[PARTidMat1;PARTidMat2];

%save PARTidMiss data
copyfile(STUDY1.PARTidMissed_filename, STUDY3.PARTidMissed_filename);
fid = fopen(STUDY3.PARTidMissed_filename,'a');
fprintf(fid,'%g, %g, %f\n', PARTidMiss2');
fclose(fid);
 
%update PARTdataMat info
if STUDY1.Gauge == STUDY2.Gauge
    PARTdataMat = [PARTdataMat1;PARTdataMat2]; %combine if same number of columns
else %otherwise pad smaller matrix with zeros
    fprintf('%s\n%s\n%s\n', 'Warning: Number of PARTdataFlds between STUDY1 and STUDY2 do not match.','Likely one study employs gauge data while other does not','PARTdataMat from study without gauge data will be padded with NaNs');
    if STUDY1.Gauge == 1 %study 1 has more columns
        PARTdataMat2 = [PARTdataMat2 nan(size(PARTdataMat2,1),size(PARTdataMat1,2)-size(PARTdataMat2,2))];
    else %study 2 has more columns
        PARTdataMat1 = [PARTdataMat1 nan(size(PARTdataMat1,1),size(PARTdataMat2,2)-size(PARTdataMat1,2))];
    end
    PARTdataMat = [PARTdataMat1;PARTdataMat2];
end

%save PARTdataMiss data
if STUDY1.Gauge == STUDY2.Gauge
    %note this might only work for windows,  need to generalize or add more
    %options for other systems
    copyEval = sprintf('copy %s+%s %s /B', STUDY1.PARTdataMissed_filename, STUDY2.PARTdataMissed_filename, STUDY3.PARTdataMissed_filename);
    system(copyEval);
else
    if STUDY1.Gauge == 1 %more columns in study1 than study2
        copyfile(STUDY1.PARTdataMissed_filename, STUDY3.PARTdataMissed_filename);
        fid = fopen(STUDY2.PARTdataMissed_filename); %add STUDY 1 PARTdataMiss
        PARTdataMiss2 = textscan(fid,'%f32 %f32 %f32','Delimiter',','); %read in STUDY 2 PARTdataMiss and pad it
        fclose(fid);
        PARTdataMiss2 = [PARTdataMiss2 nan(size(PARTdataMiss2,1), size(PARTdataMat1,2)-size(PARTdataMiss2,2))];
        fid = fopen(STUDY3.PARTdataMissed_filename,'a');
        fprintf(fid,'%g, %g, %g %u %u %u %g %g %g %g\n',PARTdataMiss2');
        fclose(fid);
    else
        fid = fopen(STUDY1.PARTdataMissed_filename); %add STUDY 1 PARTdataMiss
        PARTdataMiss1 = textscan(fid,'%f32 %f32 %f32','Delimiter',','); %read in STUDY 2 PARTdataMiss and pad it
        fclose(fid);
        PARTdataMiss1 = [PARTdataMiss1 nan(size(PARTdataMiss1,1), size(PARTdataMat2,2)-size(PARTdataMiss1,2))];
        fid = fopen(STUDY3.PARTdataMissed_filename,'w');
        fprintf(fid,'%g, %g, %g %u %u %u %g %g %g %g\n',PARTdataMiss1');
        fclose(fid);
        %note this might only work for windows,  need to generalize or add more
        %options for other systems
        copyEval = sprintf('TYPE %s>>%s', STUDY2.PARTdataMissed_filename, STUDY3.PARTdataMissed_filename);
        system(copyEval);
    end
end

%% redo PEAK data
%copy PeakMat data from STUDY1 which will remain unchanged
copyfile(STUDY1.PeakMat_filename, STUDY3.PeakMat_filename);

%get peak2 data and rewrite and append to file
fip =fopen(STUDY3.PeakMat_filename,'a');
STUDY=STUDY2;
clearlastchunk_PEAKMat=true;
rows_set2read=[];
q=char(39);  %quote character
%set up command to write out the PEAK mat data, 
fwrite4command =strcat('fwrite(fip,PEAKMat',q,',',q,'single',q,');');  
%the command is exec for each chunk read into memory
command_2do=strcat(...
    'update_idx = cell(1,length(inst1_set)); for ii= 1:size(New_PartIDs(:,1),1); update_idx{ii}=PEAKMat(:,PEAKFlds.INSTID)==New_PartIDs(ii,1); end; for ii= 1:size(New_PartIDs(:,1),1);  PEAKMat(update_idx{ii},PEAKFlds.PARTID)=PEAKMat(update_idx{ii},PEAKFlds.PARTID)+New_PartIDs(ii,3); PEAKMat(update_idx{ii},PEAKFlds.INSTID)= New_PartIDs(ii,2); end;',...
    fwrite4command);
peak_commhelper_script;
fclose(fip);

%set Peak Size info for new study info data structure
STUDY3.NumPkRows=STUDY1.NumPkRows+STUDY2.NumPkRows;
STUDY3.NumPkCols=STUDY1.NumPkCols;

PEAK = [PEAK1 PEAK2];

%% redo INST data
INST = INST1;
for i=1:length(INST2)
    findinst2 = inst2_set == INST2(i).InstID;
    if any(inst2_set == INST2(i).InstID)
        findinst1 = inst1_set(findinst2); 
        findinstrow = find([INST.InstID] == findinst1);
        peakpartid = [PEAK.PARTID];
        peakinst = [PEAK.INSTID];
        maxpartid = max(peakpartid(peakinst == findinst1));
        INST(findinstrow).LastPartID = maxpartid;
    else        
        INST = [INST INST2(i)];
    end
end;


%% Now save merged study info
STUDY=STUDY3;
fprintf('INFO, merge study,saving study info to %s\n',STUDY.Name);

save(STUDY.Name,'STUDY');
save(STUDY.DataFile,'INST','PARTidMat','PARTdataMat','PARTidFlds','PARTdataFlds','PEAK','PEAKFlds','-v7.3');
%return


function merge_study_inst(append_ornot)
%
% Call as merge_study_inst
%
% PFR 2015-09 
%   a simple merge, change inst ids in study2 and just concatenate data
%   Note, this code assumes that the INST ids in study1 are numbered in
%   sequence, so that the first INST id in study 2 will be renumbered in
%   the new merged study with a number
%   starting from the last INST id in study1
%   if append_ornot = 0 then user knows that all instids in both studies are
%   different (and all experiments are unique) and can simply append data without rewriting anything
%   if append_ornot =1 then there are overlapping instids (but not actually overlapping experiments) and new instids
%   need to be made

global YAADA STUDY DATADEF  runbatch  STUDY1 STUDY2 STUDY3 PARTdataMat PARTidMat PARTdataMat1 PARTidMat1 PARTdataMat2 PARTidMat2 PARTidMiss2 INST INST1 INST2 INST3 PEAK1 PEAK2 PEAK3 PEAK PARTidFlds PARTdataFlds1 PARTdataFlds; 

%% change inst data
if nargin<1 append_ornot=0;end;
if (append_ornot)
    instid_2start=0;
else
    instid_2start=INST1(end).InstID;  %get starting inst id for study 2
end;


for i=1:length(INST2)
        INST2(i).InstID =INST2(i).InstID+instid_2start; 
end;

%% change PART data
PARTidMat2(:,PARTidFlds.INSTID) = PARTidMat2(:,PARTidFlds.INSTID)+instid_2start;
PARTidMiss2(:,PARTidFlds.INSTID) = PARTidMiss2(:,PARTidFlds.INSTID)+instid_2start;
INST=[INST1 INST2];

%if not redoing inst need to redo partid2 so don't have identical particles
% lastPART1 = max([PARTidMat1(end,PARTidFlds.PARTID) PARTidMiss1(end,PARTidFlds.PARTID)]); %find last PARTid1
% if instid_2start == 0
% PARTidMat2(:,PARTidFlds.PARTID) = PARTidMat2(:,PARTidFlds.PARTID)+lastPART1;
% PARTidMiss2(:,PARTidFlds.PARTID) = PARTidMiss2(:,PARTidFlds.PARTID)+lastPART1;
% end

PARTidMat= [PARTidMat1;PARTidMat2];

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
        PARTdataMiss2 = cell2mat(textscan(fid,'%f32 %f32 %f32','Delimiter',',')); %read in STUDY 2 PARTdataMiss and pad it
        fclose(fid);
        PARTdataMiss2 = [PARTdataMiss2 nan(size(PARTdataMiss2,1), size(PARTdataMat1,2)-size(PARTdataMiss2,2))];
        fid = fopen(STUDY3.PARTdataMissed_filename,'a');
        fprintf(fid,'%g, %g, %g %u %u %u %g %g %g %g\n',PARTdataMiss2');
        fclose(fid);
    else
        fid = fopen(STUDY1.PARTdataMissed_filename); %add STUDY 1 PARTdataMiss
        PARTdataMiss1 = cell2mat(textscan(fid,'%f32 %f32 %f32','Delimiter',',')); %read in STUDY 2 PARTdataMiss and pad it
        fclose(fid);
        PARTdataMiss1 = [PARTdataMiss1 nan(size(PARTdataMiss1,1), size(PARTdataMat2,2)-size(PARTdataMiss1,2))];
        fid = fopen(STUDY3.PARTdataMissed_filename,'w');
        fprintf(fid,'%g, %g, %g %u %u %u %g %g %g %g\n',PARTdataMiss1');
        fclose(fid);
        copyEval = sprintf('TYPE %s>>%s', STUDY2.PARTdataMissed_filename, STUDY3.PARTdataMissed_filename);
        system(copyEval);
    end
end

%% change peak data
%merge peak row index and peak mats
%open('MergeTestBin','w');

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
% %NOTE we write it transpose so it will append ok 
command_2do=strcat('PEAKMat(:,PEAKFlds.INSTID)=PEAKMat(:,PEAKFlds.INSTID)+instid_2start;',...
                    fwrite4command);
% command_2do=strcat('PEAKMat(:,PEAKFlds.INSTID)=PEAKMat(:,PEAKFlds.INSTID)+instid_2start;',...
%     'PEAKMat(:,PEAKFlds.PARTID)=PEAKMat(:,PEAKFlds.PARTID)+lastPART1;',fwrite4command);
peak_commhelper_script;
fclose(fip);

%set Peak Size info
STUDY3.NumPkRows=STUDY1.NumPkRows+STUDY2.NumPkRows;
STUDY3.NumPkCols=STUDY1.NumPkCols;

%update PEAK structure info (the skeleton info about PEAK mat 
for i=1:length(PEAK2)
    PEAK2(i).rowcnt=PEAK2(i).rowcnt+PEAK1(end).rowcnt;
    PEAK2(i).INSTID=PEAK2(i).INSTID + instid_2start;   
%     PEAK2(i).PARTID = PEAK2(i).PARTID + lastPART1;
end;

PEAK = [PEAK1 PEAK2];

%% Now save merged study info
STUDY=STUDY3;
fprintf('INFO, merge study,saving study info to %s\n',STUDY.Name);
save(STUDY.Name,'STUDY');

if STUDY1.Gauge == 1
    PARTdataFlds = PARTdataFlds1;
end

save(STUDY.DataFile,'INST','PARTidMat','PARTdataMat','PARTidFlds','PARTdataFlds','PEAK','PEAKFlds','-v7.3');

%return


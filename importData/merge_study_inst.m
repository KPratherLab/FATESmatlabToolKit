function merge_study_inst(append_ornot)
%
% Call as merge_study_inst
% merge_study_load must be run before this
% PFR 2015-09 
%   a simple merge, change inst ids in study2 and just concatenate data
%   if append_ornot = 0 then user knows that all instids in both studies are
%   different (and all experiments and instids are unique) and can simply append data without rewriting anything
%   if append_ornot =1 then there are overlapping instids (but not actually overlapping experiments) and new instids
%   need to be made
% 
% Updated for final version
% added flexibilty if data structures don't match
% Camille Sultana 2016

global FATES STUDY DATADEF  runbatch  STUDY1 STUDY2 STUDY3 PARTdataMat PARTdataFlds PARTidMat PARTidFlds PARTidMat2 INST INST1 INST2 PEAK1 PEAK2 PEAK PEAKFlds PARTidFlds2 num1 num2 PARTidMissed PARTmisseddataFlds; 

%% change inst data
if nargin<1 append_ornot=0;end;
if ~(append_ornot)
    instid_2start=0;
else
    instid_2start=max(INST1.InstID);  %get starting inst id for study 2
end;

for i=1:length(INST2)
        INST2(i).InstID =INST2(i).InstID+instid_2start; 
end;
INST=[INST1 INST2];
STUDY3.LastInstID = max(INST.InstID);
clear global INST1 INST2

%% change PART data
%update id data for study2 hit PART data
PARTidMat((num1+1):end,PARTidFlds.INSTID) = PARTidMat((num1+1):end,PARTidFlds.INSTID)+instid_2start;
clear global PARTidMat2 PARTidMat1 PARTidFlds2 PARTidFlds1 num1 num2

%update PARTidMiss data
PARTidMissed(:,PARTidFlds.INSTID) = PARTidMissed(:,PARTidFlds.INSTID)+instid_2start;
fdataid = fopen(STUDY3.PARTidMissed_filename,'a+');
fwrite(fdataid,PARTidMissed','double');
fclose(fdataid);
clear global PARTidMissed PARTidFlds2

%% change peak data
%update INSTID in STUDY2 peakmat exernal file

%get peak2 data and rewrite and append to file
fip =fopen(STUDY3.PeakMat_filename,'a');
STUDY=STUDY2;
PEAK = PEAK2;
clearlastchunk_PEAKMat=true;
rows_set2read=[];
q=char(39);  %quote character
%set up command to write out the PEAK mat data, 
fwrite4command =strcat('fwrite(fip,PEAKMat',q,',',q,'single',q,');');  
% %NOTE we write it transpose so it will append ok 
command_2do=strcat('PEAKMat(:,PEAKFlds.INSTID)=PEAKMat(:,PEAKFlds.INSTID)+instid_2start;',...
                    fwrite4command);
peak_commhelper_script;
fclose(fip);

%set Peak Size info
STUDY3.NumPkRows=STUDY1.NumPkRows+STUDY2.NumPkRows;
STUDY3.NumPkCols=STUDY1.NumPkCols;

%update PEAK structure info (the skeleton info about PEAK mat 
for i=1:length(PEAK2)
    PEAK2(i).rowcnt=PEAK2(i).rowcnt+PEAK1(end).rowcnt;
    PEAK2(i).INSTID=PEAK2(i).INSTID + instid_2start;   
end;

PEAK = [PEAK1 PEAK2];

clear global PEAK1 PEAK2 PEAKFlds1 PEAKFlds2

%% Now save merged study info
STUDY=STUDY3;
fprintf('INFO, merge study,saving study info to %s\n',STUDY.Name);
save(STUDY.NameFull,'STUDY');

save(STUDY.DataFile,'INST','PARTidMat','PARTdataMat','PARTidFlds','PARTdataFlds','PEAK','PEAKFlds','PARTmisseddataFlds','-v7.3');

clear global STUDY3 STUDY1 STUDY2 
%return


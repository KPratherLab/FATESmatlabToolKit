function open_study
% OPEN_STUDY opens a FATES dataset
% Call as OPEN_STUDY,  
% Prequisit, call init_study and make_study to create the 'STUDY' data structure with needed
%   file information
% Open_study opens the study stored in STUDY.DataFile

global STUDY 


%PFR added info msg
fprintf('INFO, open study: about to work on %s: \n ',STUDY.Name)

if ~ischar(STUDY.NameFull) || isempty(STUDY.NameFull)
  error('Expecting string for path to Study file');
end

load(STUDY.NameFull);

if (isfield(STUDY,'DataFile'))
  fprintf('INFO, open study, about to load %s: \n ',STUDY.DataFile)
  if (exist(STUDY.DataFile,'file'))
        load(STUDY.DataFile);
        fprintf('INfO, open study, data loaded');
  else
    fprintf('ERROR, open study, no data file found to open %s \n',STUDY.DataFile);
    fprintf('   check information in STUDY structure, or perhaps you need to run make study\n');
    error(' terminating');
  end;
else
    fprintf('ERROR, open study, no data file name in STUDY \n');
    fprintf('   check pk2 directory in STUDY structure: %s\n',STUDY.PK2Dir);
    fprintf('   or perhaps you need to run make study \n');
    error(' terminating');
end;


return;


function open_study
% OPEN_STUDY opens a YAADA dataset
% Call as OPEN_STUDY,  
% Prequisit, call init_study and make_study to create the 'STUDY' data structure with needed
%   file information

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2002 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

%PFR  2014, all previous versions are not relevant as now all the data 
%  and information are stored in matlab structures and need to just be
%  loaded (without chunks)

% Fizzah Bajwa, Jonathan O. Allen  2002-06-07

% Adapted from Prakash V. Bahave's STARTUP2

% Added warning if user attempts to open YAADA v1 database
% JOA  2005-10-06

% Moved path initialization to startup.m
% JOA  2008-02-05

% Commented lines where global variables are cleared
% Alberto Cazorla  2010-09-17
% PFR 2015-4-01  streamlined version of YAADA, 

global YAADA STUDY 

%PFR  
StudyName = STUDY.Name;

%PFR added info msg
fprintf('INFO, open study: about to work on %s: \n ',StudyName)

if ~ischar(StudyName) || isempty(StudyName)
  error('Expecting string for StudyName');
end

load(StudyName);

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


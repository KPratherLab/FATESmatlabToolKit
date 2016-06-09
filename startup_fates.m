 % STARTUP starts up Matlab workspace for FATES progam and data processing
% session. Because it is a script, and because certain variables are
% defined to be global, it will create the important data structures used
% throughout FATES.  
%
% From a matlab prompt, run as 
%  >startup_fates  
%
% The program will ask the user for certain file path information.
% 
% It can also be run in a batch script as, for example:  
%    matlab -nojvm -nosplash -nodisplay -r 
%         "runbatch=1;startup_yaada;init_study('studyname');make_ydb;quit"
%
% Note that in batch, you may also need to set IN the code below;
%            runbatch=1   
%and if the startup_info.mat file is not present
%            yaada_folder = '/home/userid/yaada';
%            yaada_data= '/home/userid/mydata';
%  see also make_ydb for default settings related to batch run
%
% Note you can also change the number of peak rows to read in at atime by
%  entering a number when prompted.
%  1,000,000 is a reasonable default part(ie chunk) size for 1G-5G machines
%  In batch, which can be useful for testing or optimizing PEAK access,
%  try:
% >    startup_yaada;  
% >   YYAADA.NUMROWS_TOREADIN=1000000;
% But in batch this won't change the startup_info mmat file
%
% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  13 Oct 00
% moved bulk of code to opendb
% JOA  2002-06-07

% Restructured to use startup.mat
% Moved program parameters to YAADA in startup.mat
% Moved database parameters to STUDY in [study].mat (see initdb)
% JOA  2008-02-05

% Version changed to 2.10a (not official)
% Some modification in the folder structure. Recognizes mac, unix and pc
% PC version asks for PERL
% Alberto Cazorla 2010-09-17

% PFR 2015-4-01  streamlined version of YAADA, 
%  extensive reworking of YAADA; replaced yaada table types with structures
%      for instrument data, or just matrices for part and peak data
%    These are called, INST, PARTMat, PEAKMat; 
%    PEAKMat only holds part of all the PEAK data, which is now in a binary
%      file called PEAK_*the-study-name*.bin (if PEAK data is not large
%      PEAKMat will be able to hold all the data, depending on YAADA.NUMROWS_TOREADIN
%      value (see below).  PEAKMat is accessed with helper scripts (see
%      peak_commhelper_script and documentation)
%    PEAK is now a structure that only holds 'snapshots' of the PEAKMat
%    data, it is not used, but kept just in case
%    Also,     removed id-class code (peakid, specid,etc) to just use doubles
%              added informational messages throughout initialization
%              process. Data is not stored in 'chunks' and there are 
%              no more matlab class objects are used.  If PEAK matrix file 
%              is too large for memory, it is read in by 'parts' or
%              'chunks', but those parts are just sets of rows out of the
%              one file
%  CMS 2016 more streamlining of YAADA (removed references to SPEC
%  matrices/structures)
%    PARTdataMat and PARTidMat now only hold data for "hit" particles to reduce demands on
%    memory.  Data for "missed" particles is held in a text file called
%    PART*Missed_*the-study-name*.txt
%    Particle data is split between two matrices to save memory.  Particle
%    time and id information is held in a double precision matrix PARTidMat, while
%    all other particle data is held in a single precision matrix. 

clear all

global YAADA STUDY DATADEF INST runbatch PEAK PEAKMat PEAKFlds PARTidMat PARTidFlds PARTdataMat PARTdataFlds


runbatch=0; %default
%Note, in the case runbatch is 1, the yaada folder and data needs to be
%correct
% yaada_folder='c:\yaada_v210_prv3';
% yaada_data  =yaada_folder;


%Now if the startup info exists already use it
if exist('startup_yaada_info.mat','file')
  fprintf('INFO, startup_yaada, loading a pre-existing startup_yaada_info.mat file \n');
  load startup_yaada_info;

else  %if startup info not exist, need to ask for it
  % defaults
  YAADA.Version = '1.0 PR';

  %PFR the following fields are useful to enter the default paths for the
  %yaada_folder, where yaada code sits, and yaada_data, where raw data sits
  %NOTE, everything uses forward slash '/' b/c the '\' might be
  %interpreted as a control character when manipulating strings
  filesep2use='/';   %NOTE windows and unix will both work with /
                    %but the input command in matlab treats \ as an escape
                    %so I use / everywhere 
  if exist('yaada_folder')
     yaada_folder=strrep(yaada_folder,'\',filesep2use);
  else yaada_folder='';  %default is empty string
  end;
  %PFR, should data be associated with study? in initstudy script?
  if exist('yaada_data')
     yaada_data  =strrep(yaada_data,'\',filesep2use);
  else yaada_data='';  %default is empty string
  end;

  if isempty(yaada_folder)
     if (runbatch==0)
        yaada_folder_entered = input(sprintf('Enter the high level path for yaada folders: [%s] ',yaada_folder),'s');
        yaada_folder         = strrep(yaada_folder_entered,'\',filesep2use);
        %fprintf('You entered yaada folder: %s  \n',yaada_folder)
        if (strfind(yaada_folder_entered,'\'))
            fprintf(' note that all "\\" will be changed to "/" \n');
        end;
     end
     if isempty(yaada_folder)
         error('ERROR: The yaada folder is empty \n');
     end;
     if ~exist(yaada_folder);
          error('ERROR: The yaada folder does not exist - exiting \n');
     end;
  end
  if isempty(yaada_data)
      if (runbatch==0)
         yaada_data_entered = input(sprintf('Enter the high Level directory where your data is: [%s] ',yaada_data),'s');
         yaada_data          = strrep(yaada_data_entered,'\',filesep2use);
         %fprintf('You entered data folder: %s  \n',yaada_data)
      end
      %note, either runbatch is 1, and we cant ask the user for it,
      % or runbatch is 0,user was asked for it, but nothing was entered
     if isempty(yaada_data)
         error('ERROR: The data folder is empty \n');
     end;
     if ~exist(yaada_data);
          error('ERROR: The data folder does not exist - exiting \n');
     end;
  end

  %set up directories
  YAADA.yaada_folder = yaada_folder;  
  YAADA.yaada_data = yaada_data;  
  YAADA.MainDir    = yaada_folder;
  YAADA.TempDir    = strcat(yaada_folder,filesep2use,'tmp');
  YAADA.UserDir    = yaada_folder;  %your own scripts
  
  %set up 
  YAADA.MaxMZ = 350;
  YAADA.DeltaMZ = 0.5;
  YAADA.SpecGrav = 1.3;
%   YAADA.Verbose = 1;
%   YAADA.MaxMSView = 1000;

%PFR added this if  logic for batch executn, otherwise defaults are used
% note, to avoid defaults, you coud enter start up info so that its saved
% in a startup_info mat file
if (runbatch==0)
% 
  % YAADA Dir
  s = strrep(input(sprintf('Main YAADA directory?     [%s] ',YAADA.MainDir),'s'),'\',filesep2use);
  if ~isempty(s)
    YAADA.MainDir = s;
  end
  if ~exist(YAADA.MainDir,'dir');
    if ~mkdir(YAADA.MainDir);
      error('Cannot create directory %s',YAADA.MainDir);
    end
  end

  % User Dir
  s = strrep(input(sprintf('User program directory?   [%s] ',YAADA.UserDir),'s'),'\',filesep2use);
  if ~isempty(s)
    YAADA.UserDir = s;
  end
  if ~exist(YAADA.UserDir,'dir');
    if ~mkdir(YAADA.UserDir);
      error('Cannot create directory %s',YAADA.UserDir);
    end
  end

  % Temp Dir
  if ~exist(YAADA.TempDir,'dir');
    if ~mkdir(YAADA.TempDir);
      error('Cannot create directory %s',YAADA.TempDir);
    end
  end

  %PFR, the peak data is 1 big file and is now read in by chunks of rows
  % The more rows are read in at a time the better for performance of
  % YAADA, but it might slow down the PC
  % For a PC w/<5GB memory 1 to 5M rows should be OK,
  YAADA.NUM_PEAKROWS_TOREADIN= 1000000;  
  s = input(sprintf('Number of PEAK matrix rows to read in at a time [%d] ',...
                    YAADA.NUM_PEAKROWS_TOREADIN));
  if ~isempty(s)
      YAADA.NUM_PEAKROWS_TOREADIN = s  %s should be read as double
  end
  
  save('startup_yaada_info','YAADA');
  clear s;
end %of runbatch check


end  %checking for startup and getting folder with yaada

  
% add packages and user programs to search path
addpath(genpath(YAADA.MainDir));
addpath(genpath(YAADA.UserDir));

% Prints YAADA version and copyright notice
% Do not modify the next lines, license.txt, or credits.txt.

%PFR next lines only informational,
%type(fullfile(YAADA.MainDir,'license.txt'));
%type(fullfile(YAADA.MainDir,'credits.txt'));

return


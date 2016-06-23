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
%         "runbatch=1;startup_fates;init_study('studyname');make_ydb;quit"
%
% Note that in batch, you may also need to set IN the code below;
%            runbatch=1   
%and if the startup_info.mat file is not present
%            fates_folder = '/home/userid/fates';
%            fates_data= '/home/userid/mydata';
%  see also make_ydb for default settings related to batch run
%
% Note you can also change the number of peak rows to read in at atime by
%  entering a number when prompted.
%  1,000,000 is a reasonable default part(ie chunk) size for 1G-5G machines
%  In batch, which can be useful for testing or optimizing PEAK access,
%  try:
% >    startup_fates;  
% >   YFATES.NUMROWS_TOREADIN=1000000;
% But in batch this won't change the startup_info mmat file
%
% FATES - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  13 Oct 00
% moved bulk of code to opendb
% JOA  2002-06-07

% Restructured to use startup.mat
% Moved program parameters to FATES in startup.mat
% Moved database parameters to STUDY in [study].mat (see initdb)
% JOA  2008-02-05

% Version changed to 2.10a (not official)
% Some modification in the folder structure. Recognizes mac, unix and pc
% PC version asks for PERL
% Alberto Cazorla 2010-09-17

% PFR 2015-4-01  streamlined version of FATES, 
%  extensive reworking of FATES; replaced fates table types with structures
%      for instrument data, or just matrices for part and peak data
%    These are called, INST, PARTMat, PEAKMat; 
%    PEAKMat only holds part of all the PEAK data, which is now in a binary
%      file called PEAK_*the-study-name*.bin (if PEAK data is not large
%      PEAKMat will be able to hold all the data, depending on FATES.NUMROWS_TOREADIN
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
%  CMS 2016 more streamlining of FATES (removed references to SPEC
%  matrices/structures)
%    PARTdataMat and PARTidMat now only hold data for "hit" particles to reduce demands on
%    memory.  Data for "missed" particles is held in a text file called
%    PART*Missed_*the-study-name*.txt
%    Particle data is split between two matrices to save memory.  Particle
%    time and id information is held in a double precision matrix PARTidMat, while
%    all other particle data is held in a single precision matrix. 

clear all

global FATES STUDY DATADEF INST runbatch PEAK PEAKMat PEAKFlds PARTidMat PARTidFlds PARTdataMat PARTdataFlds


runbatch=0; %default
%Note, in the case runbatch is 1, the fates folder and data needs to be
%correct
% fates_folder='c:\fates_v210_prv3';
% fates_data  =fates_folder;

filesep2use='/';   %NOTE windows and unix will both work with /
%but the input command in matlab treats \ as an escape
%so I use / everywhere

%Now if the startup info exists already use it
if exist('startup_fates_info.mat','file')
    fprintf('INFO, startup_fates, loading a pre-existing startup_fates_info.mat file \n');
    load startup_fates_info;
    
else  %if startup info not exist, need to ask for it
    % set up default structure
    FATES.Version = '1.0 CMS';
    FATES.MaxMZ = 350;
    FATES.DeltaMZ = 0.5;
    FATES.SpecGrav = 1.3;
    
    %PFR added this if  logic for batch executn, otherwise defaults are used
    % note, to avoid defaults, you coud enter start up info so that its saved
    % in a startup_info mat file
    if (runbatch==0)
        
        % FATES Dir
        s = strrep(input('Main FATES directory?  ','s'),'\',filesep2use);
        if ~isempty(s)
            FATES.MainDir = s;
        end
        if ~exist(FATES.MainDir,'dir');
            if ~mkdir(FATES.MainDir);
                error('Cannot create directory %s',FATES.MainDir);
            end
        end
        
        % User Dir
        s = strrep(input('User program directory?  ','s'),'\',filesep2use);
        if ~isempty(s)
            FATES.UserDir = s;
        end
        if ~exist(FATES.UserDir,'dir');
            if ~mkdir(FATES.UserDir);
                error('Cannot create directory %s',FATES.UserDir);
            end
        end
        
        %Data Dir
        fates_data_entered = input('Enter the high Level directory where your data is:  ','s');
        fates_data          = strrep(fates_data_entered,'\',filesep2use);
        if isempty(fates_data)
            error('ERROR: The data folder is empty \n');
        end;
        if ~exist(fates_data);
            error('ERROR: The data folder does not exist - exiting \n');
        end;
        FATES.fates_data = fates_data;
        
        % the peak data is 1 big file and is read in by chunks of rows
        % The more rows are read in at a time the better for performance of
        % FATES, but it might slow down the PC
        % For a PC w/<5GB memory 1 to 5M rows should be OK,
        FATES.NUM_PEAKROWS_TOREADIN= 1000000;
        s = input(sprintf('Number of PEAK matrix rows to read in at a time [%d] ',...
            FATES.NUM_PEAKROWS_TOREADIN));
        if ~isempty(s)
            FATES.NUM_PEAKROWS_TOREADIN = s;  %s should be read as double
        end
        
        save('startup_fates_info','FATES');
        clear s fates_data fates_data_entered filesep2use;
    end %of runbatch check
    
    
end  %checking for startup and getting folder with fates

  
% add packages and user programs to search path
addpath(genpath(FATES.MainDir));
addpath(genpath(FATES.UserDir));

% Prints FATES version and copyright notice
% Do not modify the next lines, license.txt, or credits.txt.


return


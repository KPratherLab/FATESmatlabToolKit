function digest_tw04(PklDir,Pk2Dir)
% DIGEST_TW04 creates PK2 files from TasWare 2004 data files
% Call as DIGEST_TW04(PklDir,Pk2Dir)
% where PklDir is the directory of .pkl, .sem, .sef, .set .dam and .inst data files
%       Pk2Dir is the directory to write .pk2 files
% All files in PklDir and its subdirectories will be processed.
%
% Data directories can contain 1 instrument file for the entire directory, 
% or 1 instrument file for each .pkl file.  If there is one instrument 
% file in a directory, the data are copied to every .pk2 and the base file 
% name can be anything.  If there one instrument file for each .pkl, the 
% base file name must match the .pkl file.
% 
% This function digests raw data files created by the 2000 version of
% UCR data acquision and analysis software by Tas Dienes using the
% perl script tw00.pl.

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2002 Arizona State University

% Jonathan O. Allen  2003-05-14

% Added LaserPower field
% JOA  2003-05-14

% Added .set and light scattering functionality
% Ryan Moffet 2004-06-01

% This function was removed from the version 2.1
% It's in place again to allow TW04 digesting
% It calls tw04.pl
% PERL is called in a different way for PCs or unix/mac computers
% Alberto Cazorla 2010-09-17

% Function was rewritten for speed
% Camille Sultana 2014-04-27

global YAADA 

if ~isword(PklDir)
    error('Expecting word for PklDir');
end
if ~isword(Pk2Dir)
    error('Expecting word for Pk2Dir');
end

PerlDir = fullfile(YAADA.MainDir,'import');
DirList = sort_struct(dir(PklDir));
% handle file entry
if length(DirList) == 1
    PklDir = '';
end

%find subdirectories
isub = [DirList(:).isdir]; %returns logical vector whether name is a subfolder
nameFolds = {DirList(isub).name}; %get list of all subfolders
for i = 1:length(nameFolds)
        if strcmp(nameFolds(i),'.') || strcmp(nameFolds(i),'..');
            % ignore 
        else
            % process files in subdirectories
            digest_tw04(fullfile(PklDir,nameFolds{i}),Pk2Dir);
        end
end

%find pkl files
PklList = dir(fullfile(PklDir,'*.pkl')); %get list of all pkl files

if isempty(PklList)
    %done
else
    InstList = dir(fullfile(PklDir,'*.inst')); %find inst file
    if isempty(InstList)
        fprintf('Cannot digest %s; missing .inst file\n',PklDir);
        return;
    end
    InstFile = fullfile(PklDir,InstList(1).name);
    SemList = dir(fullfile(PklDir,'*.sem')); %find sem file
    if isempty(SemList)
        fprintf('Cannot digest %s; missing .sem file\n',PklDir);
        return;
    end
    SemFile = fullfile(PklDir,SemList(1).name);
    SetList = dir(fullfile(PklDir,'*.set')); %find set file
    if isempty(SetList)
        fprintf('Cannot digest %s; missing .set file\n',PklDir);
        return;
    end   
    SetFile = fullfile(PklDir,SetList(1).name);
    SefList = dir(fullfile(PklDir,'*.set')); %find sef file (optional)
        if isempty(SefList) 
            % sef is optional
            SefFile = '0';
        else
            SefFile = fullfile(PklDir,SefList(1).name);
        end
    if length(PklList) == 1
        PklFile = fullfile(PklDir,PklList(1).name);
        if ispc
            eval(sprintf('! %s %s %s %s %s %s %s %s',fullfile(YAADA.PerlExeDir,'perl'),fullfile(PerlDir,'tw04.pl'),InstFile,PklFile,SemFile,SetFile,SefFile,Pk2Dir));
        else
            eval(sprintf('! perl %s %s %s %s %s %s %s',fullfile(PerlDir,'tw04.pl'),InstFile,PklFile,SemFile,SetFile,SefFile,Pk2Dir));
        end
    else
        for i = 1:length(PklList)
            PklFile = fullfile(PklDir,PklList(i).name);
            if ispc
                eval(sprintf('! %s %s %s %s %s %s %s %s',fullfile(YAADA.PerlExeDir,'perl'),fullfile(PerlDir,'tw04.pl'),InstFile,PklFile,SemFile,SetFile,SefFile,Pk2Dir));
            else
                eval(sprintf('! perl %s %s %s %s %s %s %s',fullfile(PerlDir,'tw04.pl'),InstFile,PklFile,SemFile,SetFile,SefFile,Pk2Dir));
            end
        end
    end
end

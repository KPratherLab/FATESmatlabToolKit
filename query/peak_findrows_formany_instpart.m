function rows_set2use = peak_findrows_formany_instpart(instid_set,partid_set)
% a 'helper' function to read PEAK data file 
% This will find the set of rows in the PEAK Matrix file for a 
% given list of inst ids and particle ids  and return the RANGE of rows that
% would contain this set of particles.  This function is intended to speed
% up retrieval of data from the external PEAK file, by preventing having
% to read in the entire file when this is not necessary.  This works because 
% it is expected the set of particles are more or less continuous, or not spread out unusually.  
% For example,if you list 2 inst/particle ids, one inst/particle
% that occurs early in PEAK table and an inst/particle that is later, the range
% will cover all rows in between).
%  
%  Input: INSTID_2use is a Nx1 array of INSTids (ie 1 to many inst ids)
%         partid_set  is a Nx1 array of particle id values.  
%         Together each instid-partid pair (row) makes up the unique
%         identifier for a particle
%
% Output: the first and last rows from the PEAK file structure so that all data
% for the input ids is contained in between
% NOTE:  This function only finds approximately where the relevant rows are
% in PEAK Matrix file, b/c the PEAK structure only contains a rough index
% information into the matrx file.  (It will include all the relevant rows
% and more. You will not lose data) This allows you to to readin a smaller set of rows
% from the PEAK structure file and should be alot better than reading the whole
% file  (However, for small files it may not be necessary, and for
% small files that fit in memory, it's really not necessary to run this)
% 
% NOTE2: this function iterates over inst ids and calls peak_findrows_set for
% each unique inst id value found
%
%  PFR 7/2015 created
%

global PEAK PEAKMat PEAKFlds STUDY 
%partsize=min(STUDY.NumPkRows,100000);  %1Mrows x 9 double ~ 72Mb

if (nargin<2 || (length(instid_set)~=length(partid_set)) )
    disp('Call peak_findrows_formany_instpart(instid_set,partid_set)');
    error('  ..where instid set and partid set are equal size arrays');
end;

instid_uq = unique(instid_set); %find unique values in instid_set

curr_nr =size(PEAKMat,1);
if (curr_nr>=STUDY.NumPkRows)  
       %then current in memory PEAK Matrix is everything
       %just do find on peakmat
     disp('All peak data for study currently loaded in PEAK matrix. No need to load data from external file.');
     rows_set2use = [];
     return
else
    temprows=[];
    for i=1:length(instid_uq)
        part2use= partid_set(instid_set == instid_uq(i)); %get partids that match each unique instid
        if (iscell(part2use)) part2use=part2use{1}; end %get array not cell
        temprows=[temprows; peak_findrows_set(instid_uq(i),[min(part2use) max(part2use)])]; %find rows in PEAK file that contain data for these particles
    end;
    rows_set2use=[min(temprows(:,1)) max(temprows(:,2))]; %find rows in PEAK file that contain data for all input particles
end;

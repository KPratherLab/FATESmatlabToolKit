function rows_set2use = peak_findrows_set(INSTID_2use,part_range)
% This will find the set of rows in the PEAK Matrix file for a 
% given range of particle ids and a single inst id
%  
% Input: INSTID_2use is a single InstID (a single integer)
%        part_range is a 2 element array of particle id values
%        eg  [1000 2000] will find the rows in PEAK Matrix file that
%        contain particles with PID >=1000 and <=2000 with INSTID =
%        INSTID_2use.
% Output: an index to the first and last rows from the PEAK structure 
% so that all the peak data for the input id values is contained between
% the indexed rows
% NOTE:  This function only finds approximately where the relevant rows are
% in PEAK Matrix file, b/c the PEAK structure only contains a rough index
% information into the matrx file.  (It will include all the relevant rows
% and more. You will not lose data) This allows you to to readin a smaller set of rows
% from the PEAK structure file and should be alot better than reading the whole
% file  (However, for small files it may not be necessary, and for
% small files that fit in memory, it's really not necessary to run this)
%
%  PFR 2/2015 created
%

global PEAK PEAKMat PEAKFlds STUDY 
%partsize=min(STUDY.NumPkRows,100000);  %1Mrows x 9 double ~ 72Mb
  
curr_nr =size(PEAKMat,1); %get number of rows
if (curr_nr==STUDY.NumPkRows)  
       %then current in memory PEAK Matrix is everything
       %just do find on peakmat
     disp('All peak data for study currently loaded in PEAK matrix. No need to load data from external file.');
else
    lowerbound = INSTID_2use == [PEAK.INSTID] & part_range(1) >= [PEAK.PARTID];
    lowerbound = find(lowerbound,1,'last');
    upperbound = INSTID_2use == [PEAK.INSTID] & part_range(2) <= [PEAK.PARTID];
    upperbound = find(upperbound,1);
    chunk_set2use=[lowerbound upperbound];
end;


%todo need to considre instid=x but part > < something
%consider 4 args, instid,relop, and return blkset then get the actuall peak
%data
if (0)
    find_results=[];
    fid         =fopen(STUDY.PeakMat_filename,'r');
    
    
    peakmatraw =fread(fid,numcols*partsize,'single');
    pkloadcnt= pkloadcnt+1;
    numpkrows  =size(peakmatraw,1)/numcols;
    if (numpkrows > 0)
        PEAKMat    =reshape(peakmatraw,numcols,numpkrows);
        clear peakmatraw;
        PEAKMat = PEAKMat';
        fprintf('INFO, loaded PEAKMat, %i with  %i rows \n',pkloadcnt,numpkrows);
        %now peak is mattrix and runquery 2 should use PEAKMat
        run_find();
    end;
    
end;
%now convert chunks to rows in file to read
rows_set2use=[PEAK(chunk_set2use(1)).rowcnt PEAK(chunk_set2use(2)).rowcnt];

end

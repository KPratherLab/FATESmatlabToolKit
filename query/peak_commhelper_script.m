% -----------------------------------------------------------------
% This is a helper script to process the PEAK table in parts (or chunks).
% Note that 'parts' are just sets of rows to read in at a time,
%  as given in the YAADA.NUM_PEAKROWS_TOREADIN variables 
%    (or if not present, the default is  1,000,000)
%
% PFR created 4/2015
% To use this script call it from matlab, or from your own script
%
%
% Example:  this will read in 1000 rows and leave it in PEAKMat in memory
%    rows_set2read =[1 1000];
%    command_2do='';
%    clearlastchunk_PEAKMat=false;
%    peak_commhelper_script
%
% Example: this will read in the whole data, 
%   but if there's more rows in data file than
%   the num-peak-rows-2readin from YAADA structure, it will only keep the
%   last set readin
%    rows_set2read   =[];
%    command_2do     ='';
%    clearlastchunk_PEAKMat=false;
%    peak_commhelper_script
%
% Example: this will read in the whole data and find indices where part
% id=8 and leave PEAKMat empty, but save results in myPkmat
%  (note that your command can be test by first getting a few rows from 
%    PEAKMat (see above example), then running it in matlab, then lastly running 
%    peak_commhelper_script)
%
%    rows_set2read =[];
%    myPkmat       =[];
%    command_2do   ='i=find(PEAKMat(:,2)==8);
%                    myPkmat=[myPkmat;PEAKMat(i,:)];';
%    clearlastchunk_PEAKMat=true;
%    peak_commhelper_script
%
% In Detail:
% 1. Your calling session, or calling scripts, must have a variable called 'command_2do'
%    This variable holds a text string, this string is the command that
%        is run on each part of the peak data matrix
%    If this variable does not exist then this program will readin data and
%    do nothing (so you may still see messages)
%
%   Note the text string can be a set of matlab statements or the name of 
%   a script containing matlab statements.  Those statements process
%   PEAKMat, and this helper scripts loads in parts of the data from the 
%   binary file into the PEAKMat matrix in memory, before calling the script. 
%   The command_2do can be empty and 
%   the first part of the data matrix will be read into memory
%
% 2. Also you must have a variable called 
% rows_set2read, an array of 2 numbers:[first-row last-row] to readin  
%which is used to filter peak rows to run on only a part of
% whole peak matrix file, this can be empty and it will run on all rows
%
% 3. Optional 
%    clearlastchunk_PEAKMat is a logical variable (true or false)
%    to indicate if PEAKMat should be 
%    cleared out (set to empty matrix) before finishing this script.
%    This may be useful to save memory, or it can be unused and this helper
%    script will leave PEAKMat with the last chunk read in
% -----------------------------------------------------------------


%function [find_results,find_results_index]= peak_findhelper(find_command2do,rows_set2read)
%  A helper function to search the PEAK Matrix file
%  This will perform file reads. load parts of the file,
%  and perform the find command on those parts
%  Parameter: find_string is  the find command to execute
%  Returns  : concatenated results of find, which is index set
%
% PFR 2/2015  created

global PEAK PEAKMat PEAKFlds STUDY YAADA
comm_results      =[];
pkloadcnt         =0;    %count number of chunks read in
pk_bytesperelem   =4;    %assume PEAK Matrix file was written out as single

%guess a good partition size to read in at a time
%need to check if peak is already in memory
if (isfield(YAADA,'NUM_PEAKROWS_TOREADIN'))
    NUM_ROWS_TOREADIN = YAADA.NUM_PEAKROWS_TOREADIN;
else
    NUM_ROWS_TOREADIN = 1000000;
end;
partsize=min(STUDY.NumPkRows,NUM_ROWS_TOREADIN);  %1Mrows x 9 double ~ 72Mb

%if (nargin==1)
if (~exist('rows_set2read') || length(rows_set2read)==0)
%       fprintf('INFO, peak comm helper script, setting rows_set2read to all rows\n');
       rows_set2read=[1 PEAK(end).rowcnt];  %default is to read whole file
elseif (~(length(rows_set2read)==2))
         error('ERROR, the row set should be 2 numbers equal to first/last rows');
end; 

%check the curr number of row in PEAKMat
%if the current set of rows in memory,ie in PEAKMat Matrix, 
% is everything from file,then just evaluate command
curr_nr =size(PEAKMat,1);
if (curr_nr>=STUDY.NumPkRows)  
       fprintf('INFO, peak comm helper, PEAKMat all %i rows in memory already\n',curr_nr);
       
%        %evaluate command on all data or just the rows set requested
%        if (rows_set2read(1)==1 && rows_set2read(2) == PEAK(end).rowcnt)
%           % run_cmd();  %run find on whole peakmat
%           eval(command_2do);
%           fprintf('INFO, peak comm helper script, results, %i th chunk \n',1);
%                     
%        else   %run the find commnd on the rows set range given
%           tempMat   =PEAKMat;   %save current PEAKMat to avoid any side effects
%           PEAKMat   =PEAKMat(rows_set2read(1):rows_set2read(2),:);
%           %run_cmd(); 
%           eval(command_2do);
%           fprintf('INFO, peak comm helper script,rows range filtered, results, %i th chunk\n',1);
%           PEAKMat   =tempMat;  %restore PeakMat
%        end;
else
     fprintf('INFO, peak comm helper, will read rows from %i to %i\n',rows_set2read(1),rows_set2read(2));
     %set up loop through 'chunks' of PEAKMat that are in matrix file
     fid         = fopen(STUDY.PeakMat_filename,'r');
     numcols     = STUDY.NumPkCols;  %or STUDY.NumPkCols
     pkloadcnt   = 0;
     
     %get 1st row to read and byte offset into file
     %NOTE, 4 bytes assumes the peakmat file was written as singles in
     %parse_part when imported
     offset_2start    = (rows_set2read(1)-1)*numcols*pk_bytesperelem;  %9 elems per row,4bytes per elem     
     last_rowread     = rows_set2read(1)-1;  %keep track of last matrix file row read in
     fseek(fid,offset_2start,-1);    %subtract 1 from beg of file, goto offset byte
     
     while (~feof(fid) && last_rowread < rows_set2read(2))  
          %%total elements to read in is partsize*numcols
          peakmatraw   = fread(fid,numcols*partsize,'single');  
          pkloadcnt    = pkloadcnt+1;
          numpkrows    =size(peakmatraw,1)/numcols;  %calculate numrows read in
          last_rowread = last_rowread+numpkrows;     %update row counts
                
       if (numpkrows>0)  %make sure something was read in
          %Note, the data is read into a Nx1 vector, the first 9 elements 
          % are the first row, and so on, so we must reshape and transpose
          % (note also that reshape works by filling in first column
          %   of the LHS matrix result, 
          %   so the reshape has to first use numcols as the 1st index)
          PEAKMat    =reshape(peakmatraw,numcols,numpkrows);
          %get memory usage
          if ispc
          mem_stats  = memory;
          end;
          clear peakmatraw;
          PEAKMat    = PEAKMat';
          if (last_rowread>rows_set2read(2))  %if we read past rqsted range
                row_overcnt =last_rowread - rows_set2read(2);
                PEAKMat     =PEAKMat(1:end-row_overcnt,:);
                last_rowread=rows_set2read(2);
          end;
          %fprintf('INFO, loaded PEAKMat chunk %i, with  %i file rows, %i mat rows, last row is %i \n',...
          %          pkloadcnt,numpkrows,size(PEAKMat,1),last_rowread);
          if ispc
               fprintf('%f Mb, %i \n', mem_stats.MemUsedMATLAB/1000000, last_rowread);
          end
          eval(command_2do);  

%          fprintf('INFO, peak helper script, processed, %i th chunk\n',pkloadcnt);
%          comm_results      =[comm_results;temp_results];
       end;  %end if numpkrws>0
     
     end; %end while
     fprintf('INFO, pkcomm helper, Summary: loaded PEAKMat in %i chunks, with  %i total file rows \n',...
            pkloadcnt,numpkrows);
     if (exist('clearlastchunk_PEAKMat') && clearlastchunk_PEAKMat==true)
          PEAKMat=[];
     end
     fprintf('INFO, pkcomm helper, left PEAKMat with %i mat rows, last row read in %i \n',...
             size(PEAKMat,1),last_rowread);

        
end


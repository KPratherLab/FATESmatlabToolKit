function [missedData, hitData] = read_part(HitPartFile,MissedPartFile)
% READ_PART reads in the files containing the particle data for the hit (currently .set) 
% and missed (currently .sem) particles into matrices. These matrices are used
% by parse_part to populate the PART data tables.  This file may be altered to  
% accomodate other file formats (other than .sem and .set).  Efforts should
% be made when altering the code to use the fastest data import MATLAB
% function appropiate, as this step of study creation can often be a
% bottleneck. 

global missedFlds hitFlds

%% get part data for missed particles sem files
if ~isempty(MissedPartFile)
    fim = fopen(MissedPartFile);
    if fim ~= -1 %change code from here to else (line 18) statement if using a different file format
        missedData = textscan(fim, '%*f %*f %f %f %s %f %f %f %f','TreatAsEmpty',{'1.#J'},'Delimiter',','); %note first two columns are ignored, last seven columns are read in 
        fclose(fim);
        %change string date to num
        missedData{:,missedFlds.TIME} = (datenum(missedData{:,missedFlds.TIME}));
        missedData = cell2mat(missedData); %make matrix
    else
        missedData = []; %get size
    end
else
    missedData = [];
end

%% get part data for hit particles set files
if ~isempty(HitPartFile)
    fit = fopen(HitPartFile);
    if fit ~= -1 %change code from here to else (line 34) statement if using a different file format
        hitData = textscan(fit, '%*f %*f %f %f %s %f %f %f %f','TreatAsEmpty',{'1.#J'},'Delimiter',','); %note first two columns are ignored, last seven columns are read in 
        fclose(fit);
        %change date to string
        hitData{:,hitFlds.TIME} = (datenum(hitData{:,hitFlds.TIME}));
        hitData = cell2mat(hitData); %make matrix
    else
        hitData = []; %do not change
    end
else
    hitData = []; %do not change
end
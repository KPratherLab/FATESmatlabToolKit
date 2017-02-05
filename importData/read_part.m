function [missedData, hitData] = read_part(HitPartFile,MissedPartFile)
% READ_PART reads in the files containing the particle data for the hit (currently .set) 
% and missed (currently .sem) particles into matrices (missedData and hitData). These matrices are used
% by parse_part to populate the PART data tables.  This file may be altered to  
% accomodate other file formats (other than .sem and .set).  Efforts should
% be made when altering the code to use the fastest data import MATLAB
% function appropiate, as this step of study creation can often be a
% bottleneck. 

global missedFlds hitFlds

%% CODE FOR PRATHER ATOFMS FILE FORMAT (.sem and .set)
%comment out or change if not using Prather ATOFMS file format
% 
%get part data for missed particles sem files
if ~isempty(MissedPartFile)
    fim = fopen(MissedPartFile);
    if fim ~= -1 %change code from here to fclose if using a different file format
        missedData = textscan(fim, '%*f %*f %f %f %s %f %f %f %f','TreatAsEmpty',{'1.#J'},'Delimiter',','); %note first two columns are ignored, last seven columns are read in 
        fclose(fim);
        %don't change anything below here
        %change string date to num
        missedData{:,missedFlds.TIME} = (datenum(missedData{:,missedFlds.TIME},'dd-mmm-yyyy HH:MM:SS'));
        missedData = cell2mat(missedData); %make matrix
    else
        missedData = []; %do not change
    end
else
    missedData = []; %do not change
end

%get part data for hit particles set files
if ~isempty(HitPartFile)
    fit = fopen(HitPartFile);
    if fit ~= -1 %change code from here to fclose if using a different file format
        hitData = textscan(fit, '%*f %*f %f %f %s %f %f %f %f','TreatAsEmpty',{'1.#J'},'Delimiter',','); %note first two columns are ignored, last seven columns are read in 
        fclose(fit);
        %don't change anything below here
        %change date to string
        hitData{:,hitFlds.TIME} = (datenum(hitData{:,hitFlds.TIME},'dd-mmm-yyyy HH:MM:SS'));
        hitData = cell2mat(hitData); %make matrix
    else
        hitData = []; %do not change
    end
else
    hitData = []; %do not change
end

%% CODE FOR ALABAMA FILE FORMAT (HitParticles.txt and MissedParticles.txt)
%comment out or change if not using ALABAMA file format
% 
% %get part data for missed particles DetectedParticles.txt
% if ~isempty(MissedPartFile)
%     fim = fopen(MissedPartFile);
%     if fim ~= -1 %change code from here to fclose if using a different file format
%         missedData = textscan(fim, '%s %f','Delimiter',{'\t'}); %reads in data from all detected particle file
%         fclose(fim);
%         missedData{:,missedFlds.TIME} = (datenum(missedData{:,missedFlds.TIME},'dd.mm.yyyy HH:MM:SS')); %change string date to num
%         missedData = cell2mat(missedData); %make matrix
%     else
%         missedData = []; %do not change
%     end
% else
%     missedData = []; %do not change
% end
% 
% %get part data for hit particles HitParticles.txt
% if ~isempty(HitPartFile)
%     fit = fopen(HitPartFile);
%     if fit ~= -1 %change code from here to fclose if using a different file format
%         hitData = textscan(fit,'%s %f %*[^\n]','Delimiter','\t'); %read in data from hit particles (particle and spectra info)
%         fclose(fit);
%         hitData{:,hitFlds.TIME} = (datenum(hitData{:,hitFlds.TIME},'dd.mm.yyyy HH:MM:SS')); %change date string to num
%         hitData = cell2mat(hitData); %make matrix
%     else
%         hitData = []; %do not change
%     end
% else
%     hitData = []; %do not change
% end
% 
% %eliminate hit particle data from missed particle matrix
% if ~isempty(missedData) && ~isempty(hitData)
%     missedData = setdiff(missedData, hitData, 'rows');
% end

%% CODE FOR COMMERCIAL ATOFMS FILE FORMAT (.sem and .set)
%comment out or change if not using commerical ATOFMS file format
% % 
% % get part data for missed particles sem files
% if ~isempty(MissedPartFile)
%     fim = fopen(MissedPartFile);
%     if fim ~= -1 %change code from here to fclose if using a different file format
%         missedData = textscan(fim, '%*f %f %*f %s','Delimiter',','); %note first two columns are ignored, last seven columns are read in 
%         fclose(fim);
%         %don't change anything below here
%         %change string date to num
%         missedData{:,missedFlds.TIME} = (datenum(missedData{:,missedFlds.TIME},'mm/dd/yyyy HH:MM:SS'));
%         missedData = cell2mat(missedData); %make matrix
%     else
%         missedData = []; %do not change
%     end
% else
%     missedData = []; %do not change
% end
% 
% %get part data for hit particles set files
% if ~isempty(HitPartFile)
%     fit = fopen(HitPartFile);
%     if fit ~= -1 %change code from here to fclose if using a different file format
%         hitData = textscan(fit, '%*f %*s %f %*f %f %s','Delimiter',','); %note first two columns are ignored, last seven columns are read in 
%         fclose(fit);
%         %don't change anything below here
%         %change date to string
%         hitData{:,hitFlds.TIME} = (datenum(hitData{:,hitFlds.TIME},'mm/dd/yyyy HH:MM:SS'));
%         hitData = cell2mat(hitData); %make matrix
%     else
%         hitData = []; %do not change
%     end
% else
%     hitData = []; %do not change
% end
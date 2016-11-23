function [PARTidMissed, PARTdataMissed] = load_missed(PidMissedFile, PdataMissedFile)
%LOAD_MISSED loads in missed particle data
%Call as [PARTidMissed PARTdataMissed] = load_missed(PidMissedFile, PdataMissedFile)
%INPUTS
%PidMissedFile: fullfile path to the binary file containing the pids and time of all
%missed particles in the study
%PdataMissedFile: fullfile path the binary file containing the particle
%data for all missed particles in the study
%These variables take up a lot of memory, which is why it is not held in the 
%study at all the times.
%OUTPUT
%PARTidMissed: particle id matrix following same format as PARTidMat, each
%row corresponds to data for a single particle
%PARTdataMissed: particle data matrix following same format as PARTdataMat,
%each row corresponds to data for a single particle

global PARTidMat PARTdataMat STUDY
PARTidMissed = [];
PARTdataMissed = [];

%open up PidMissedFile
tmp = fopen(PidMissedFile);
if tmp ~= -1
    PARTidMissed = fread(tmp,'double');
    fclose(tmp);
    
    %shape matrix
    numRows = length(PARTidMissed)/size(PARTidMat,2);
    PARTidMissed = reshape(PARTidMissed, [size(PARTidMat,2) numRows]);
    PARTidMissed = PARTidMissed';
end


%open up PdataMissedFile
tmp = fopen(PdataMissedFile);
if tmp ~= -1
    PARTdataMissed = fread(tmp,'single');
    fclose(tmp);
    
    %shape matrix
    numRows = length(PARTdataMissed)/size(PARTdataMat,2);
    numCols = STUDY.NumPartCols;
    PARTdataMissed = reshape(PARTdataMissed, [numCols numRows]);
    PARTdataMissed = PARTdataMissed';
end

end
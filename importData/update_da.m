function update_da
% UPDATE_DA updates Da in all PART chunks
% Call as UPDATE_DA

% YAADA - Software Toolkit to Analyze Single-Particle Mass Spectral Data
%
% Copyright (C) 1999-2000 California Institute of Technology
% Copyright (C) 2001-2008 Arizona State University
% Copyright (C) 2008 Jonathan O. Allen

% Jonathan O. Allen  28 Jun 00

% Revised for global STUDY
% JOA  2008-02-06
% Revised for streamlined Yaada 
% PFR/CMS 2016

global INST PARTidFlds PARTidMat PARTdataMat PARTdataFlds

%PFR
uniq_iids=unique([INST(:).InstID]); %find unique InstIDs
for i =1:length(uniq_iids)
    InstIdx=find([INST(:).InstID]==uniq_iids(i),1); %determine row for InstID
    fprintf('INFO, update da, instidx: %i',InstIdx)
    
    func2call=char(INST(InstIdx).DaCalibFunction); %retrieve function to calculate da
    
    if ~isempty(func2call)
        %evaluate function using calibration parameters from inst table and
        %particle velocities
        func_result=feval(func2call,INST(InstIdx).DaCalibParam,PARTdataMat(PARTidMat(:,PARTidFlds.INSTID)==uniq_iids(i),PARTdataFlds.VELOCITY));
        
        %add Da as a field if not already there
        if ( ~ismember('DA',fieldnames(PARTdataFlds)) )
            lastfld=length(fieldnames(PARTdataFlds));
            PARTdataFlds.DA=lastfld+1;
        end;
        PARTdataMat(PARTidMat(:,PARTidFlds.INSTID)==uniq_iids(i),PARTdataFlds.DA)=func_result;   %add DA to PartMat table
        
    else fprintf('WARN, update da, no da calib function to call');
    end;
end;

return
function loadDATA_alabama
% Call as loadDATA_alabama
% loadDATA_alabama identifies data files to be loaded into the FATES study
% and then writes them to the FATES study variables.
% loadDATA_alabama is called by make_study

global STUDY nameDET nameHIT PARTdataMat procDATE INST formatHit MZValue Polarity nameINST

%% locate data files
findData_alabama(STUDY.RawDir);
save(STUDY.DataList,'nameDET','nameHIT','procDATE','nameINST');

%% set up some values used in parse_part_alabama
%format string to read in HitParticle file
formatHit = '%s %f';
for i = 1:500
    formatHit = strcat(formatHit,' %f32');
end

%mz values for a spectra
MZValue = single([-1*(1:250) 1:250]);
MZValue = MZValue';
numMZ = length(MZValue);

%spectra polarity values
Polarity = single([zeros(numMZ/2,1); ones(numMZ/2,1)]);

%% parse and write particle and inst and peak data
single(PARTdataMat);
for i = 1:length(nameHIT)
    [InstID,NewInst] = parse_inst(nameINST{i});
    disp(nameDET{i})
    parse_part_alabama(nameHIT{i},nameDET{i},InstID);
end

clearvars -global formatHit MZValue Polarity nameHIT nameDET hitData missedData numFldsPARTdata peakFldsNAME spectraNAME
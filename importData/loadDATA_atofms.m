function loadDATA_atofms
% Call as loadDATA_atofms
% loadDATA_atofms identifies data files to be loaded into the FATES study
% and then writes them to the FATES study variables.
% loadDATA_atofms is called by make_study

global STUDY nameSET namePKL nameSEM nameINST PARTdataMat procDATE

%locate data files
findData(STUDY.RawDir);
save(STUDY.DataList,'nameSEM','nameSET','namePKL','nameINST','procDATE');

%parse and write particle and inst and peak data
single(PARTdataMat);
for i = 1:length(nameSET)
    [InstID,NewInst] = parse_inst(nameINST{i});
    disp(nameSET{i})
    parse_part(nameSET{i},nameSEM{i},namePKL{i},InstID);
end
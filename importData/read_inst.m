function instText = read_inst(InstFile)
% read_inst reads the file containing instrument/experiemental conditions into the variable instText
% read_inst is called by parse_inst which uses instText to fill data into
% the FATES structure INST. 
% Call as instText = read_inst(InstFile)
% The input InstFile is a cell containing the full path to the .inst file
% Currently read_inst is written to digest a text file. 
% The output instText should be a Nx2 cell array.  The first column contains
% the fieldname (variable name) that will be saved to the INST structure
% and the second column contains the matching variable data.  instText must maintain
% this format to remain compatible with parse_inst. 

fid = fopen(InstFile{1},'rt'); %open up .inst file
instText = textscan(fid,'%s','Delimiter','\n'); %read in inst file
instText = instText{1}(~cellfun('isempty',instText{1})); %get rid of empty lines
instText = strtrim(instText); %get rid of trailing and leading white spaces
instText = instText(~strncmp(instText,'%',1)); %get rid of lines that are comments
instText = regexp(instText,'\=','split'); %split string
instText = strtrim(instText); %get rid of any trailing/leading white spaces
fclose(fid);

end
function [NegResponse,PosResponse] = get_NONint_spectrum_SUM(PID,MZbins,ResponseType,Polarity)
% GET_NONINT_SPECTRUM_SUM returns spectrum matrix summed over the bins provided 
% Call as [NegResponse,PosResponse] = get_NONint_spectrum_SUM(PID,MZbins,ResponseType,Polarity)
% Where PID is set of paritcle identifiers stored as a nx2 matrix
%   PID(:,1) = InstID
%   PID(:,2) = PartID
% MZbins are the bin EDGES to be used in histcounts (not center of bin).  The default MZbins is
%  .  This gives gives half integer bins centered over 0.5 and 1 
% values ie 0.5 (0.25-0.75), 1 (0.75-1.25), 1.5 (1.25-1.75)
% The ResponseType can be any column in PEAK. 'Area' is the default ResponseType.
% Polarity specifies the spectrum polarity as
%   Polarity = 0 - negative spectra
%   Polarity = 1 - positive spectra
%   Polarity = 2 - negative and positive spectra (default)
%
% NegResponse and PosResponse are matrices with columns matching hit particles
% and rows represent summed response over each bin. Remember when plotting NegResponse
% or PosResponse against mz use the mid point of each bin, not the input
% MZbins which are the bin edges

global FATES

%% check inputs
%check setup
if nargin < 1 || nargin > 4
  error('Call as [NegResponse,PosResponse] = get_NONint_spectrum_SUM(PID,MZbins,ResponseType,Polarity)');
end

%check PID is correct size (Nx2 matrix)
if (~size(PID,2)==2)  
  error('Invalid PID');
end

%check MaxMZ
if exist('MZbins','var')
  if ~isvector(MZbins) || length(MZbins) < 2
    error('Expecting vector for MZbins');
  end
else
  MZbins = 0.25:0.5:300.5; %this gives half integer bins centered over 0.5 and 1 values ie 0.5 (0.25-0.75), 1 (0.75-1.25), 1.5 (1.25-1.75)
end

%check ResponseType
if ~exist('ResponseType','var')
  ResponseType = 'RelArea';
else
  if ~ischar(ResponseType);
    error('Expecting word for ResponseType');
  end
end

%check Polarity
if ~exist('Polarity','var')
  Polarity = 2;
else
  if Polarity < 0 || Polarity > 2
    error('Expecting 0 (negative), 1 (positive), or 2 (both) for Polarity');
  end
end
if Polarity == 2 && nargout ~= 2
  error('For dual polarity, call with 2 output arguments');
end

%% get spectra and bin

%get raw spectra
Spectrum = get_spectrum(PID,Polarity,ResponseType);

%mem_stats=memory; for cking memory on windows 
fprintf('INFO, getintspectrum sum, matrix size: %i X %i \n',size(Spectrum));

%set up output variables
NumPart = size(PID,1);
switch Polarity
  case 0 
    NegResponse = single(zeros(size(MZbins,2)-1,NumPart));
    PosResponse = [];
  case 1
    NegResponse = [];
    PosResponse = single(zeros(size(MZbins,2)-1,NumPart));
  case 2
    NegResponse = single(zeros(size(MZbins,2)-1,NumPart));
    PosResponse = single(zeros(size(MZbins,2)-1,NumPart));
end

%get integer spectra for each particle
for i = 1:NumPart
    if ~isempty(Spectrum{i,1}) %make sure particle has spectrum data
        MZ = Spectrum{i,1}(:,1); 
        Response = Spectrum{i,1}(:,2);
        
        %neg spectra
        if ~isempty(NegResponse)
            idx = MZ(:)<-MZbins(2) & MZ(:)>= -MZbins(end); %find neg mz
            allNEG = Response(idx); %get responses for neg mz
            [~,~,mzIDX] = histcounts(MZ(idx),fliplr(-MZbins));
            if any(mzIDX)
                mzIDX = length(MZbins)-mzIDX;
                NegTMP = accumarray(mzIDX,allNEG); %add all peaks within an integer unit together
                NegResponse(1:length(NegTMP),i) = NegTMP;
            end
        end
        
        % positive spectra
        if ~isempty(PosResponse)
            idx = MZ(:)<MZbins(end) & MZ(:)>= MZbins(2); %find pos mz
            allPOS = Response(idx); %get responses for pos mz
            [~,~,mzIDX] = histcounts(MZ(idx),MZbins);
            if any(idx)
                PosTMP = accumarray(mzIDX,allPOS); %add all peaks within an integer unit together
                PosResponse(1:length(PosTMP),i) = PosTMP;
            end
        end
    end
end

if nargout == 1
    % return only requested polarity
    if Polarity == 1
        NegResponse = PosResponse;
    end
    
end

return
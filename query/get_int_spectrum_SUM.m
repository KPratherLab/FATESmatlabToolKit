function [NegResponse,PosResponse] = get_int_spectrum_SUM(PID,MaxMZ,ResponseType,Polarity)
% GET_INT_SPECTRUM_SUM returns spectrum matrix summed over integral m/z 
% Call as [NegResponse,PosResponse] = get_int_spectrum_SUM(PID,MaxMZ,ResponseType,Polarity)
% Where PID is set of paritcle identifiers stored as a nx2 matrix
%   PID(:,1) = InstID
%   PID(:,2) = PartID
% MaxMZ is the upper limit of m/z range; The default MaxMZ is FATES.MaxMZ.
% The ResponseType can be any column in PEAK. 'Area' is the default ResponseType.
% Polarity specifies the spectrum polarity as
%   Polarity = 0 - negative spectra
%   Polarity = 1 - positive spectra
%   Polarity = 2 - negative and positive spectra (default)
%
% NegResponse and PosResponse are matrices with columns matching hit particles
% and rows represent integral MZ.  PosResponse(23,1) is the sum of areas of peaks with MZ =
% 22.5-23.5 for the first hit particle.
% The PosResponse columns span m/z = 1 to MaxMZ, the NegResponse columns span m/z = -1 to -MaxMZ. 

global FATES

%check setup
if nargin < 1 || nargin > 4
  error('Call as [NegResponse,PosResponse] = get_int_spectrum_SUM(PID,MaxMZ,ResponseType,Polarity)');
end

%check PID is correct size (Nx2 matrix)
if (~size(PID,2)==2)  
  error('Invalid PID');
end

%check MaxMZ
if exist('MaxMZ','var')
  if ~isnumeric(MaxMZ) || length(MaxMZ) > 1
    error('Expecting scalar for MaxMZ');
  end
else
  MaxMZ = FATES.MaxMZ;
end

%check ResponseType
if ~exist('ResponseType','var')
  ResponseType = 'RelArea';
else
  if ~ischar(ResponseType);
    error('Expecting string for ResponseType');
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

%get raw spectra
Spectrum = get_spectrum(PID,Polarity,ResponseType);

%mem_stats=memory; for cking memory on windows 
fprintf('INFO, getintspectrum sum, matrix size: %i X %i \n',size(Spectrum));

%set up variables
NumPart = size(PID,1);
switch Polarity
  case 0 
    NegResponse = single(zeros(MaxMZ,NumPart));
    PosResponse = [];
  case 1
    NegResponse = [];
    PosResponse = single(zeros(MaxMZ,NumPart));
  case 2
    NegResponse = single(zeros(MaxMZ,NumPart));
    PosResponse = single(zeros(MaxMZ,NumPart));
end

%get integer spectra for each particle
for i = 1:NumPart
    if ~isempty(Spectrum{i,1}) %make sure particle has spectrum data
        MZ = Spectrum{i,1}; 
        Response = Spectrum{i,2};
        
        %neg spectra
        if ~isempty(NegResponse)
            idx = (MZ(:)<-0.5 & MZ(:)>= -(MaxMZ+0.5)); %find neg mz
            allNEG = Response(idx); %get responses for neg mz
            if any(idx)
                NegIntMZ = ceil(abs(MZ(idx))*(1-10*eps)-0.5); %get integer mz
                NegResponse(1:max(NegIntMZ),i) = accumarray(NegIntMZ,allNEG); %add all peaks within an integer unit together
            end
        end
        % positive spectra
        if ~isempty(PosResponse)
            idx = (MZ(:)<(MaxMZ+0.5) & MZ(:)>= 0.5); %find pos mz
            allPOS = Response(idx); %get responses for pos mz
            if any(idx)
                PosIntMZ = ceil(MZ(idx)*(1+10*eps)-0.5); %get integer mz
                PosResponse(1:max(PosIntMZ),i) = accumarray(PosIntMZ,allPOS); %add all peaks within an integer unit together
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



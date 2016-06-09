function [A,B] = massCalibration(X,Y,plotFlag)
% massCalibration performs a calibration given two arrays of (1)
% points from a ATOFMS spectrum and (2) the corresponding M/Z values for
% those points. Adapted from WDR_DLL script (T Rebotier) 3/3/2016 using the
% following: The model is Y = A.(X-B)^2 where Y is the M/Z and X is the time (point number)
%We therefore seek the minimum of Q = SUM_i( (Y_i-A.(X_i-B)^2)^2 )
%for simplicity in the following formulas I omit the index i but keep in mind 
%that all sums are over different X_i's and Y_i's so only expressions in A or B purely
%can be factored out of the sums; we note Ui = (Xi-B)^2:
%Q = SUM( (Y-A.(X-B)^2)^2 )							(1)
%dQ/dA = 2 SUM ( (AU-Y)U )							(2)
%dQ/dA = 0 gives us A = SUM(YU)/SUM(U^2)				(3)
%dQ/dB = -4.A.SUM( (X-B).(AU-Y) )					(4)
%dQ/dB = 0 gives us B = SUM(X(AU-Y))/SUM(AU-Y)		(5)
%This formula cannot be solved directly but can be numerically (it converges)
%we iterate (3) and (5) starting with A =1 and B=0, after a few runs we keep the pair
%(a long and painful first order development of B(t+1) shows it 
%converges to the optimal B in the neighborhood of that optimal B...
%
%Changed (TR 7/22) to newB = average ( X - root(Y/A))
%A = slope, B = intercept
%To find M/Z using these two parameters:
% MZ = (A*([1:numPoints]-B).^2)';

if nargin == 2
    plotFlag = 0;
end

initial_A = 2E-6; % Define constants and initial values
initial_B = 0;
iterations = 1000;
A = initial_A;
B = initial_B;

for i = 1:iterations % Loop to determine A and B
    clear temp_B
    SYU = 0;SU2 = 0;newB = 0;
    U = (X-B).^2;
    U2 = U.*U;
    SYU = sum(Y.*U);
    SU2 = sum(U.*U);
    temp_B = X-(Y./A).^(0.5);
    A = SYU/SU2; % Define new A
    newB = sum(temp_B)/length(X); % Find new B, change slowly
    B = (0.9*B) + (0.1*newB);

end

%% Optional plot to visually determine whether the fit is good
if plotFlag == 1 
    new_MZ = (A*([1:15000]-B).^2)'; % calculate new m/z
    figure;
    plot([1:15000],new_MZ,'LineWidth',2); % plot new spectrum
    hold on
    plot(X,Y,'rx','LineWidth',2); % plot points
end
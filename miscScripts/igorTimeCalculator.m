function [igorTime] = igorTimeCalculator(timeString)
matlabTime = datenum(timeString);
igorTime = (matlabTime-datenum('01-JAN-1904 00:00:00'))*(3600*24);

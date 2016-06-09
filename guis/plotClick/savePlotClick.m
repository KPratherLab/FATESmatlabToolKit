%script for saving plots for multiPlotClick

if plotSave == 1
    if plotChoose == 0 %save plots automatically
        exportToPPTX('addslide');
        exportToPPTX('addpicture', q');
        savefig(fullfile(plotFolder,Fname));
    else
        plotSaveNow = input('\nSave current figure? 0-no, 1-yes  '); %save select plots
        if plotSaveNow == 1
            savefig(fullfile(plotFolder,Fname));
            exportToPPTX('addslide');
            exportToPPTX('addpicture', q');
        end
    end
end
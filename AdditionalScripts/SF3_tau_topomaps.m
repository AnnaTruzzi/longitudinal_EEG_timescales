% Clear all the space and variables
clc
clear all 

%% Load the tau data and separe it by age; 

Path2tau   = '/Users/jr4545/Desktop/Tau_data/';  
dataname   = 'tau_exploratory_v4.csv'; 

data     = readtable([Path2tau filesep dataname]); 
data.Label(strcmp(data.Label, 'Cz')) = {'E129'};
elect    = unique(data.Label); 


%% Load any dataset and transform it to fieldtrip structure to create good quality topomap
Path2EEGlab = '/Users/jr4545/Downloads/eeglab2024.2/';
Path2Field  = '/Users/jr4545/Downloads/fieldtrip-20181231/';
datapath    = '/Users/jr4545/Desktop/';
dataname    = 'BUDDY144_6mo_Res.set'; 

addpath(genpath(Path2EEGlab))
addpath(genpath(Path2Field))

EEG = pop_loadset(dataname,datapath);

EEG = pop_select(EEG, 'channel', elect);

EEG = eeglab2fieldtrip(EEG,'preprocessing');

% remove the innecesary data of the dataset 

try
EEG = rmfield(EEG,'trialinfo'); 
end 
try
EEG = rmfield(EEG,'time'); 
end
try
EEG = rmfield(EEG, 'trial');
end
try
EEG.dimord   = 'chan_freq'; 
end
EEG.freq     = 1; 

%% Introduce the name of the topomaps you'll create and the path to save them 

path2figure = Path2tau;


%% LOOP to creat the topomapts - this loop create a basic topomap in average

powtype = {'its_exploratory'}; 
ses     = [6, 9, 16]; 

for powidx  = 1:length(powtype)

    powdata = data.t; 

    
    for sesidx = 1:3
        
        dataprov = powdata(data.ses_age == ses(sesidx));
        
        for chanidx = 1:length(EEG.label)
            electpos = find(strcmp(EEG.label{chanidx}, elect)==1);
            dataprov(chanidx) = dataprov(electpos); 
        end 
        dataprov = dataprov(1:length(EEG.label));
        
        dataplot(:,sesidx) = dataprov; 

        if sesidx == 3 
            minval = 0.01; 
            maxval = 0.30; 
            for ageidx = 1:3
                            plotname = [powtype{powidx} '_' num2str(ses(ageidx)) 'mo_tau.jpeg'];

            EEG.powspctrm(:,1) = dataplot(:,ageidx);

            cfg = [];
            cfg.marker       = 'off';
            cfg.zlim         = [minval maxval];
            cfg.ylim         = [1,1];
            cfg.comment      = 'no'; 
            cfg.interactive  = 'no'; 
            cfg.colormap     = 'pink';
            cfg.colorbar     = 'SouthOutside';
 
            figure; 
            ft_topoplotTFR(cfg, EEG)
                set(gcf,...
                'units','centimeters',...
                'position', [0,0,5,5])
                           set(gca,... 
                'FontUnits','points',...
                'DefaultaxesFontName', 'Times new Roman',...
                'DefaultaxesFontSize', 11,...
                'Color', 'white')
            print([path2figure filesep plotname],'-djpeg', '-r1200')
            close all
            end
            
        end 
        
        
    end
    end 



%% LOOP to creat the topomapts - this loop create a basic topomap in average - Validation

dataname   = 'tau_exploratory_v4.csv'; 

data     = readtable([Path2tau filesep dataname]); 
data.Label(strcmp(data.Label, 'Cz')) = {'E129'};

elect    = unique(data.Label); 

powtype = {'its_validation'}; 
ses     = [6, 9, 16]; 

for powidx  = 1:length(powtype)

    powdata = data.t; 

    
    for sesidx = 1:3
        
        dataprov = powdata(data.ses_age == ses(sesidx));
        
        for chanidx = 1:length(EEG.label)
            electpos = find(strcmp(EEG.label{chanidx}, elect)==1);
            dataprov(chanidx) = dataprov(electpos); 
        end 
        dataprov = dataprov(1:length(EEG.label));
        
        dataplot(:,sesidx) = dataprov; 

        if sesidx == 3 
            minval = .01; 
            maxval = .40; 
            for ageidx = 1:3
                            plotname = [powtype{powidx} '_' num2str(ses(ageidx)) 'mo_tau.jpeg'];

            EEG.powspctrm(:,1) = dataplot(:,ageidx);

            cfg = [];
            cfg.marker       = 'off';
            cfg.zlim         = [minval maxval];
            cfg.ylim         = [1,1];
            cfg.comment      = 'no'; 
            cfg.interactive  = 'no'; 
            cfg.colormap     = 'pink';
            cfg.colorbar     = 'SouthOutside';
 
            figure; ft_topoplotTFR(cfg, EEG);
                set(gcf,...
                'units','centimeters',...
                'position', [0,0,5,5])
                           set(gca,... 
                'FontUnits','points',...
                'DefaultaxesFontName', 'Times new Roman',...
                'DefaultaxesFontSize', 11,...
                'Color', 'white')
            print([path2figure filesep plotname],'-djpeg', '-r1200')
            close all
            end
            
        end 
        
        
    end
    end 

%% LOOP to creat the topomapts - this loop create a basic topomap in average - Adults

dataname   = 'tau_adults_v4.csv'; 

data     = readtable([Path2tau filesep dataname]); 
data.Label(strcmp(data.Label, 'Cz')) = {'E129'};

elect    = unique(data.Label); 

powtype    = {'its_adults'}; 
condition  = ["Video", "EO", "EC"]; 

for powidx  = 1:length(powtype)

    powdata = data.t; 

    
    for sesidx = 1:3
        
        dataprov = powdata(strcmp(data.event, condition{sesidx}));
        
        for chanidx = 1:length(EEG.label)
            electpos = find(strcmp(EEG.label{chanidx}, elect)==1);
            dataprov(chanidx) = dataprov(electpos); 
        end 
        dataprov = dataprov(1:length(EEG.label));
        
        dataplot(:,sesidx) = dataprov; 

        if sesidx == 3 
            minval = 0.01; 
            maxval = .40; 
            for ageidx = 1:3
                            plotname = [powtype{powidx} '_' condition{ageidx} '_tau.jpeg'];

            EEG.powspctrm(:,1) = dataplot(:,ageidx);

            cfg = [];
            cfg.marker       = 'off';
            cfg.layout       = '/Users/jr4545/Desktop/BFY_A6/egi128_layout.sfp';
            cfg.zlim         = [minval maxval];
            cfg.ylim         = [1,1];
            cfg.comment      = 'no'; 
            cfg.interactive  = 'no'; 
            cfg.colormap     = 'pink';
            cfg.colorbar     = 'SouthOutside';
      
            figure; ft_topoplotTFR(cfg, EEG);
                set(gcf,...
                'units','centimeters',...
                'position', [0,0,5,5])
                      set(gca,... 
                'FontUnits','points',...
                'DefaultaxesFontName', 'Times new Roman',...
                'DefaultaxesFontSize', 11,...
                'Color', 'white')
            print([path2figure filesep plotname],'-djpeg', '-r1200')
            close all
            end
            
        end 
        
        
    end
    end 
function [A, locs] = sweepSacf()

    pars = loadParameters();
    pars.est.dur = 40;
    pars.est.noiseOff = 5;
    pars.est.nOfIts = 1;
    pars.est.onset = 0;
    pars.est.tail  = 0;
    pars.est.f     = 1000;
    pars.est.type  = 'FMsweep';
    pars.subDelay  = 0;
    pars.subCortAff = 1;
    pars.regularise = 0;
    pars.tauSubthal = 1;
    pars.tauSACF    = 2.5;

    lagSpace = binSpace(pars);

    f0    = [900, 1200, 1500];
    delta = [linspace(-600, 600, 10)];

    for i = 1:length(f0)
        pars.est.f = f0(i);
        for j = 1:length(delta)
            tt = tic;
            locs{i, j} = zeros([10, 20]);
            fprintf('Computing f0 = %.0f, delta = %.1f .', f0(i), delta(j))
            pars.est.shift = delta(j);
            parfor k = 1:10
                [ATmp(k, :),  locsTmp(k, :)] = pyThalamic(lagSpace, pars);
            end
            A{i, j}    = ATmp; 
            locs{i, j} = locsTmp; 
            fprintf('done! Time = %.1fm\n', toc(tt) / 60);  
            save('../sacfSims.mat', 'lagSpace', 'locs', 'A')
        end
    end

    %f = 1000./mean(locs' ./ (1:length(locs)));

end



function [A, locs] = pyThalamic(lagSpace, pars)
 
    % parse filename is randomized to allow parallel computations
    chart = char(['A':'Z' 'a':'z' '0':'9']);
    parseID = ['pyparse' chart(ceil(length(chart) * rand(1, 4)))];

    save([parseID, 'In.mat'], 'pars', 'lagSpace');
 
    %python = '/home/tabs/Apps/anaconda/bin/ipython --colors=NoColor';
    python = 'python';

    system([python, ' sacf.py ', parseID]);

    pyparse = load([parseID, 'Out.mat']);
    delete([parseID, 'Out.mat']);
    
    A = mean(pyparse.A, 1);    
    [~, locs] = findpeaks(flip(A), flip(lagSpace), 'NPeaks', 20);
    
end



function lagSpace = binSpace(pars)

    lagMin = 1000 / (pars.freqInterval(2));
    lagMax = 1000 / (pars.freqInterval(1));
    lagSpace = linspace(lagMax, lagMin, pars.N)';

end



function pars = loadParameters()

    pars.freqInterval = [33, 2000]; % [minFreq, maxFreq] (Hz)
    pars.N            = 250;        % number of exc/inh populations 
    pars.cortFs       = 1000;       % sample rate in cortex (samples per sec)
    pars.subCortAff = 4;     % number of afferents in the subcortical pathway
    pars.regularise = 1;     % whether to regularise (1) or not (0) the SACF
    pars.mu0        = 75;    % expected peak activity after normalisation (Hz)
    pars.tauSubthal = -1;    % subthalamic lowpass tau (ms)
    pars.solvOnset  = 1;     % Solvability onset for the delay channels (ms)
    pars.SACFGround = 0.35;  % SACF baseline ([0,1]; -1 for average signal)
    pars.tauSACF    = -1;    % sacf tau (-1 for Wiegribe's taus; ms)

end

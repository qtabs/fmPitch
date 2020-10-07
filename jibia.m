function jibia()

    addpath('/home/tabs/Cloud/Libs/matlab/distributionPlot')
        
    % First octopus.py and sacf.m to generate the result files. This script
    % is only for plotting!

    clear; close all
    fig1();  % Sweep Pitch Shift
    fig2();  % variance and up/down asymmetry
    fig3();  % spectral model
    fig4();  % temporal model
    %fig6();  % responses examples
    fig7();  % direction selectivity
    fig8();  % tau(I, h)
    fig9();  % heatmaps predictions sweep pitch shift
    fig10(); % predictions on sweep variance and up/down asymmetry
    fig11(); % parameter space exploration
    fig13(); % Brady's stimuli
    fig14(); % train sweep pitch shift  
    fig15(); % heatmaps predictions train sweep pitch shift
    fig16(); % model predictions up/down train asymmetry 
    fig17(); % Stimulus example 
    %fig18(); % Connectivity matrices
    figS1(); % Effect of presentation order on Delta p
    figS2(); % Slope of the linear fits per subject

end


function fig1() % Sweep Pitch Shift

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    swFit = 1:10;
    load('deltaPitchDB.mat');

    f0  = avgFreqs{1};
    dSw = deltas{1};
    sw  = results;

    for i = 1:length(f0)
        FSw{i} = [];
        fSw{i} = [];
        for j = 1:length(subLab)
            FSw{i}  = [FSw{i};  sw{j}.F{i}(:, :) - f0(i)];
            fSw{i}  = [fSw{i};  sw{j}.F{i}(:, swFit) - f0(i)];
        end

        x = repmat(dSw(swFit), [length(subLab) * 2, 1]);
        X = repmat(dSw(swFit), [length(subLab) * 4, 1]);

        [p, s] = polyfit(X, fSw{i}, 1);
        fitErr = sqrt(diag(inv(s.R) * inv(s.R')) .* s.normr .^ 2 ./ s.df);
        fitResid = sum((fSw{i}(:) - polyval(p, X(:))).^2); 
        r2Sw(i)  = 1 - fitResid / ((length(fSw{i}(:)) - 1) * var(fSw{i}(:))); 
        fitSwm(i) = p(1);
        errSwm(i) = fitErr(1);
        fitSwn(i) = p(2);
        errSwn(i) = fitErr(2);

    end

    fig = figure;
    for i = 1:length(f0)
        subplot(1, length(f0), i); hold off;
        x = {'-600', '', '-333', '', '-66', '66', '', '333', '', '600'};
        distributionPlot(FSw{i}, 'showMM', 4, 'xNames', x, ...
                         'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
        title(['$\bar{f} = ' int2str(f0(i)) '$\,Hz'], 'Interpreter', 'LaTeX');
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('$\Delta p$ (Hz)', 'Interpreter', 'LaTeX')
        ylim([-600, 600])
        yticks(-500:200:500);
        zeroLine = refline(0, 0);
        zeroLine.Color = 'k';
        %fitLine = refline(fitSwm(i),fitSwn(i));
        m = fitSwm(i) * (dSw(3)- dSw(1)) / (3 - 1);
        n = fitSwm(i) * dSw(1) - m;
        fitLine = refline(m, n);
        fitLine.Color = .2 * [1 1 1];
        fitLine.LineWidth = 1;
        fitLine.LineStyle = '--';
        uistack(zeroLine, 'bottom'); 
    end

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig1.svg', '-dsvg');

end



function fig2() % variance and up/down asymmetry

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    swFit = 1:10;
    trFit = 1:6;

    load('deltaPitchDB.mat');

    f0  = avgFreqs{2};
    dSw = deltas{2};

    for sIx = 1:length(subLab)
        [~, ix] = max(strcmp(subjects, subjects{sIx}));
        sw{sIx} = results{ix};
    end 

    for i = 1:length(f0)
        FSw{i} = [];
        fSw{i} = [];
        for j = 1:length(subLab)
            FSw{i}  = [FSw{i};  sw{j}.F{i}(:, :) - f0(i)];
            fSw{i}  = [fSw{i};  sw{j}.F{i}(:, swFit) - f0(i)];
        end

        x = repmat(dSw(swFit), [length(subLab) * 2, 1]);
        X = repmat(dSw(swFit), [length(subLab) * 4, 1]);

    end

    fig = figure;
    subplot(121); hold off;
    for i = 1:length(subLab)
        for j = 1:length(f0) 
            swVar(length(f0) * (i - 1) + j, :) = std(sw{i}.F{j} - f0(j));
        end
    end
    x = {'-600', '', '', '-200', '', '', '200', '', '', '600'};
    distributionPlot(swVar, 'showMM', 4, 'xNames', x, ...
                     'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    ylabel('$std(\Delta p)$ (Hz)',  'Interpreter', 'LaTeX')
    xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
    title('variance of $\Delta p$', 'Interpreter', 'LaTeX')
    ylim([0, 260])
    yticks(0:100:200);
    swVarAbs = [swVar(:, 1:5); swVar(:, 10:-1:6)];
    swDeltaMat = repmat(linspace(600, 0, 5), [size(swVarAbs, 1), 1]);
    [swVarCorr, swVarP] = corr(swDeltaMat(:), swVarAbs(:), 'type', 'spearman');
      

    subplot(122); hold off;
    fSwAll = [fSw{1}; fSw{2}; fSw{3}];
    asymmTmp = fSwAll(:, 1:5) + fSwAll(:, 10:-1:6);
    udSw = 0.25 * (asymmTmp(1:4:end,:) + asymmTmp(2:4:end,:) + ...
                   asymmTmp(3:4:end,:) + asymmTmp(4:4:end,:));
    x = {'600', '', '333', '', '66'};
    distributionPlot(udSw, 'showMM', 4, 'xNames', x, ...
                    'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (Hz)',  'Interpreter', 'LaTeX')
    title('up-down asymmetry', 'Interpreter', 'LaTeX')
    yticks(-200:100:200);
    ylim([-200, 260])
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    uistack(zeroLine, 'bottom'); 

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2*plosCBwidths(2) 1.5];
    print(fig, 'fig2.svg', '-dsvg');

    fprintf('Avg up/down asymmetry: %.0fHz\n', sqrt(mean(udSw(:).^2)));

end



function fig3() % spectral model

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('spectSims.mat');
    
    f0 = [900, 1200, 1500];
    delta = linspace(-600, 600, 10);
    sSpace = 1:size(p.f900d600, 2);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            condName = sprintf('f%.0fd%.0f', f0(i), delta(j));
            c{i}(j, :) = mean(p.(condName));
            d{i}(j, :) = mean(expFr.(condName));
        end
    end


    fig = figure;
    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(j) = sum(sSpace .* fDist);
            sweepSigm(j) = (sum(sSpace.^2 .* fDist) - sweepPeak(j)^2);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(j) = sum(sSpace .* fDist);
        end
        subplot(1, 3, i)
        hold off
        imagesc(delta, sSpace, c{i}');
        hold on
        errorbar(delta, sweepPeak, sweepSigm, 'LineWidth', 2, 'Color', 'black');

        plot(delta, expChannel, 's', 'Color', 'black', 'MarkerSize', 10);
        xticks([-600, -333.33, -66.66, 66.66, 333.33, 600])
        xticklabels([-600, -333, -66, 66, 333, 600])
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('cochlear channel index', 'Interpreter', 'LaTeX')
        ylim([25, 65])
        title(['$\bar{f} = ', sprintf('%d', f0(i)), '\,$Hz'], 'Interpreter', 'LaTeX')
    end

    legend('spectral predictions', 'experimental results')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig3.svg', '-dsvg');

    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(i,j) = sum(sSpace .* fDist);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(i,j) = sum(sSpace .* fDist);
        end
    end

    R2 = 1-sum((expChannel-sweepPeak).^2)/sum((expChannel-mean(expChannel)).^2);

    fprintf('spectral model: R2 = %.3f\n', R2);

end



function fig4() % temporal model

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    % First run sweepSacf from octopus.py to 
    % generate spectSims.mat. This script only plots the results
    n0 = 1; n1 = 4;
    f0 = [900, 1200, 1500];

    load('sacfSims.mat');
    sweepExp = load('expResAvgSw.mat');
    delta = linspace(-600, 600, 10);

    tSpace{1} = linspace(0.8, 1.53, 50);
    tSpace{2} = linspace(0.6, 1.10, 50);
    tSpace{3} = linspace(0.5, 0.84, 50);
    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            for k = 1:length(tSpace{i})
                c{i}(j, k) = 0;
                for n = n0:n1
                    [~, ix] = min((lagSpace - n * tSpace{i}(k)).^2);
                    c{i}(j, k) = c{i}(j, k) + sum(A{i, j}(:, ix));
                end
            end
        end
        c{i}(j, k) = c{i}(j, k) / ((n1 - n0) * size(A{i, j}, 1));
    end

    fig = figure;
    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (c{i}(j, :).^2 - min(c{i}(j, :).^2)) / sum(c{i}(j, :).^2 - min(c{i}(j, :).^2));
            sweepPeak(j) = sum(tSpace{i} .* fDist);
            sweepSigm(j) = (sum(tSpace{i}.^2 .* fDist) - sweepPeak(j)^2);
        end
        subplot(1, 3, i)
        hold off
        imagesc(delta, tSpace{i}, c{i}');
        hold on
        errorbar(delta, sweepPeak, sweepSigm, 'LineWidth', 2, 'Color', 'black');

        plot(delta, 1000 ./ sweepExp.fAvg(i, :), 's', 'Color', 'black', 'MarkerSize', 10);
        xticks([-600, -333.33, -66.66, 66.66, 333.33, 600])
        xticklabels([-600, -333, -66, 66, 333, 600])
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('period (ms)', 'Interpreter', 'LaTeX')
        ylim([min(tSpace{i}), 0.1 * floor(10 * max(tSpace{i}))])
        title(['$\bar{f} = ', sprintf('%d', f0(i)), '\,$Hz'], 'Interpreter', 'LaTeX')
    end

    legend('SACF predictions', 'experimental results')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig4.svg', '-dsvg');


    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (c{i}(j, :).^2 - min(c{i}(j, :).^2)) / sum(c{i}(j, :).^2 - min(c{i}(j, :).^2));
            sweepPeak(i,j) = sum(tSpace{i} .* fDist);
            expPeriod(i,j) = 1000 ./ sweepExp.fAvg(i, j);
        end
    end

    R2 = 1-sum((expPeriod-sweepPeak).^2)/sum((expPeriod-mean(expPeriod)).^2);

    fprintf('temporal model: R2 = %.3f\n', R2);

end



function fig6() % responses examples

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('spectSims.mat');

    channels = 1:100;
    cLims = [0, 100];
    fig = figure;

    load('popresDown.mat')
    timeSpace = timeSpace - 20;
    fs = int16(1 / (timeSpace(2)-timeSpace(1))); % in ms

    subplot(6, 2, 2)
    imagesc(timeSpace, channels, hUpExcit', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, 'sweep layer $\uparrow$ (excitatory)', 'Interpreter', 'LaTeX', 'Color', 'white')
    cb = colorbar('Location', 'East');
    cb.Color = 'white';
    cb.Ticks = [0, 100];
    cb.Label.String = 'firing rate (Hz)';
    cb.Label.Interpreter = 'LaTeX';
    yticks([]);

    subplot(6, 2, 4)
    imagesc(timeSpace, channels, hDownExcit', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, 'sweep layer $\downarrow$ (excitatory)', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])

    subplot(6, 2, 6)
    imagesc(timeSpace, channels, hFrequency', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, 'spectral layer', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])

    subplot(6, 2, 8)
    imagesc(timeSpace, channels, 0.1 * peripheral', cLims)
    ylabel('cochlear channels $\longrightarrow$', 'Interpreter', 'LaTeX')
    text(-45, 85, 'auditory periphery', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])
    set(gca,'YDir','normal')

    subplot(6, 2, 10)
    plot(timeSpace((21*fs):(70*fs)), 0.001 * fStim((21*fs):(70*fs)), 'k-', 'LineWidth', 2);
    ylabel('simulus $f(t)$ (kHz)', 'Interpreter', 'LaTeX')
    xlabel('time from stimulus onset (ms)', 'Interpreter', 'LaTeX')
    xlim([-5, 80])

    subplot(6, 2, 12)
    plot(channels, sum(hFrequency) / sum(hFrequency(:)), 'w-', 'LineWidth', 2);
    xlabel('channels', 'Interpreter', 'LaTeX')


    load('popresUp.mat')
    channels = 1:100;
    cLims = [0, 100];
    timeSpace = timeSpace - 20;
    fs = int16(1 / (timeSpace(2)-timeSpace(1))); % in ms

    subplot(6, 2, 1)
    imagesc(timeSpace, channels, hUpExcit', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, '$\uparrow$ excitatory', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])

    subplot(6, 2, 3)
    imagesc(timeSpace, channels, hDownExcit', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, '$\downarrow$ excitatory', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])

    subplot(6, 2, 5)
    imagesc(timeSpace, channels, hFrequency', cLims)
    set(gca,'YDir','normal')
    text(-45, 85, 'frequency integrators', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])

    subplot(6, 2, 7)
    imagesc(timeSpace, channels, 0.1 * peripheral', cLims)
    ylabel('cochlear channels $\longrightarrow$', 'Interpreter', 'LaTeX')
    text(-45, 85, 'peripheral response', 'Interpreter', 'LaTeX', 'Color', 'white')
    yticks([])
    xlim([-5, 80])
    set(gca,'YDir','normal')

    subplot(6, 2, 9)
    plot(timeSpace((21*fs):(70*fs)), 0.001 * fStim((21*fs):(70*fs)), 'k-', 'LineWidth', 2);
    ylabel('simulus $f(t)$ (kHz)', 'Interpreter', 'LaTeX')
    xlabel('time from stimulus onset (ms)', 'Interpreter', 'LaTeX')
    xlim([-5, 80])

    subplot(6, 2, 11)
    plot(channels, sum(hFrequency) / sum(hFrequency(:)), 'w-', 'LineWidth', 2);
    xlabel('channels', 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 7];
    print(fig, 'fig6.svg', '-dsvg');

end



function fig7() % direction selectivity

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('dsi.mat')

    fig = figure;

    f0 = [900, 1200, 1500];
    delta = linspace(200/3, 600, 5);

    for i = 1:length(f0)
        for j = 1:length(delta)
            d = topDownOn.downExc(i, 11-j) - topDownOn.downExc(i, j);
            dsiDown(i,j) = d/(topDownOn.downExc(i,11-j)+topDownOn.downExc(i,j));
            d = topDownOn.upExc(i, 11-j) - topDownOn.upExc(i, j);
            dsiUp(i,j) = d/(topDownOn.upExc(i,11-j)+topDownOn.upExc(i,j));
        end
    end

    dsiUp   = fliplr(dsiUp);
    dsiDown = fliplr(dsiDown);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
    end

    subplot(121)
    hold off
    for i = 1:length(f0)
        plot(delta, dsiUp(i, :), 'Color', plosCBcolors(i, :), 'Marker', styles(i));
        hold on
    end
    ylabel('DSI$^{\uparrow}$',  'Interpreter', 'LaTeX')
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    l = legend(freqLegend, 'Interpreter', 'LaTeX', 'Location', 'SouthEast');
    %yticks(500:200:1100)
    xticks(round(delta));
    ylim([0, 1])
    xlim([0, 633])

    subplot(122)
    hold off
    for i = 1:length(f0)
        plot(delta, dsiDown(i, :), 'Color', plosCBcolors(i, :), 'Marker', styles(i));
        hold on
    end
    ylabel('DSI$^{\downarrow}$',  'Interpreter', 'LaTeX')
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    %yticks(500:200:1100)
    xticks(round(delta));
    ylim([-1, 0])
    xlim([0, 633])

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.25*plosCBwidths(2) 1.75];
    print(fig, 'fig7.svg', '-dsvg');


    for i = 1:length(f0)
        for j = 1:length(delta)
            d = topDownOff.downExc(i, 11-j) - topDownOff.downExc(i, j);
            offDsiDown(i,j) = d/(topDownOff.downExc(i,11-j)+topDownOff.downExc(i,j));
            d = topDownOff.upExc(i, 11-j) - topDownOff.upExc(i, j);
            offDsiUp(i,j) = d/(topDownOff.upExc(i,11-j)+topDownOff.upExc(i,j));
        end
    end

    offDsiUp   = fliplr(offDsiUp);
    offDsiDown = fliplr(offDsiDown);

    ratOnOff = (dsiUp(:) - dsiDown(:)) ./ (offDsiUp(:) - offDsiDown(:));
    avgRatOnOff = mean(ratOnOff);
    errRatOnOff = std(ratOnOff)/sqrt(length(ratOnOff));
    ratOnOffStr = sprintf('%.02f %s %.02f',avgRatOnOff,char(177),errRatOnOff);

    fprintf('Ratio of DSI with top-down on/off = %s\n', ratOnOffStr);

end



function fig8() % tau(I, h)

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    a = 310;
    b = 125;
    d = 0.16;

    x = 0.23 + [0.02, 0.05, 0.065, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, ...
                0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.215, 0.24, 1];
    h = linspace(0, 50, 1000)';

    y   = a * x - b;
    phi = y ./ (1 - exp(-d * y));

    dphidx  = phi .* (1 ./ y + d ./ (1 - exp(d .* y)));
    tau = 100 * dphidx ./ h;

    tau(tau > 20) = 20;

    %colmap = ([1 .31 0]' * linspace(1, 0.3, 100))';
    colmap = colormap(winter(size(tau, 2)));

    fig = figure; hold off

    for i = 1:length(x)
        plot(h, tau(:, i), 'Color', colmap(i, :), 'LineWidth', 1);
        hold on;
    end


    % Traj 1
    x = [linspace(0.23, 0.43, 40), 0.43 * ones([1, 960])];
    h = [linspace(0, 4, 300).^2, linspace(16, 40, 700)];
    y   = a * x - b;
    phi = y ./ (1 - exp(-d * y));
    dphidx  = phi .* (1 ./ y + d ./ (1 - exp(d .* y)));
    traj1 = 100 * dphidx ./ h;
    traj1(traj1 > 20) = 20;
    plot(h, traj1, 'k--', 'LineWidth', 1.5);

    % Traj 2
    x = [linspace(0.23, 0.3, 100)];
    h = [linspace(0, 4, 100)];
    y   = a * x - b;
    phi = y ./ (1 - exp(-d * y));
    dphidx  = phi .* (1 ./ y + d ./ (1 - exp(d .* y)));
    traj1 = 100 * dphidx ./ h;
    traj1(traj1 > 20) = 20;
    plot(h, traj1, 'k:', 'LineWidth', 1.5);

    % Traj 2
    x = [linspace(0.3, 0.43, 40), 0.43 * ones([1, 660])];
    h = [linspace(2, 6, 40).^2, linspace(36, 40, 660)];
    y   = a * x - b;
    phi = y ./ (1 - exp(-d * y));
    dphidx  = phi .* (1 ./ y + d ./ (1 - exp(d .* y)));
    traj1 = 100 * dphidx ./ h;
    traj1(traj1 > 20) = 20;
    plot(h, traj1, 'k-', 'LineWidth', 1.5);

    text(20, 16, 'initial state', 'Interpreter', 'LaTeX')
    text(20, 14, 'feedback modulation', 'Interpreter', 'LaTeX')
    text(20, 12, 'encoding (not modulated)', 'Interpreter', 'LaTeX')
    text(20, 10, 'encoding (pre-modulated)', 'Interpreter', 'LaTeX')
    text(20,  8, 'final state (equilibrium)', 'Interpreter', 'LaTeX')

    c = colorbar('East');
    colormap(colmap)
    c.Label.String = 'synaptic input';
    c.Label.Interpreter = 'LaTeX';
    c.Ticks = [];

    ylabel('integration time constant (ms)', 'Interpreter', 'LaTeX')
    xlabel('firing rate (Hz)', 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 3.5 2];
    print(fig, 'fig8.svg', '-dsvg');

end



function fig9() % heatmaps predictions sweep pitch shift

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('pitchPredictionsSw.mat');

    f0 = [900, 1200, 1500];
    delta = linspace(-600, 600, 10);
    sSpace = 1:size(sp.f900d600, 2);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            condName = sprintf('f%.0fd%.0f', f0(i), delta(j));
            c{i}(j, :) = mean(sp.(condName));
            d{i}(j, :) = mean(ep.(condName));
        end
    end

    fig = figure;
    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(j) = sum(sSpace .* fDist);
            sweepSigm(j) = (sum(sSpace.^2 .* fDist) - sweepPeak(j)^2);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(j) = sum(sSpace .* fDist);
        end
        subplot(1, 3, i)
        hold off
        imagesc(delta, sSpace, c{i}');
        hold on
        errorbar(delta, sweepPeak, sweepSigm, 'LineWidth', 2, 'Color', 'black');

        plot(delta, expChannel, 's', 'Color', 'black', 'MarkerSize', 10);
        xticks([-600, -333.33, -66.66, 66.66, 333.33, 600])
        xticklabels([-600, -333, -66, 66, 333, 600])
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('cochlear channel index', 'Interpreter', 'LaTeX')
        ylim([25, 65])
        title(['$\bar{f} = ', sprintf('%d', f0(i)), '\,$Hz'], 'Interpreter', 'LaTeX')
    end

    legend('FM-feedback predictions', 'experimental results')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig9.svg', '-dsvg');


    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(i,j) = sum(sSpace .* fDist);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(i,j) = sum(sSpace .* fDist);
        end
    end

    channSpan = mean(expChannel(:, end) - expChannel(:, 1));
    expAvgHz  = load('expResAvgSw.mat');
    freqSpan  = mean(expAvgHz.fAvg(:, end) - expAvgHz.fAvg(:, 1));

    R2 = 1-sum((expChannel-sweepPeak).^2)/sum((expChannel-mean(expChannel)).^2);
    err2 = sqrt(mean((expChannel(:)-sweepPeak(:)).^2));

    fprintf('FM-Feedback model (sweeps): R2 = %.3f, SME = %.3f\n', R2, err2);
    fprintf('Hz/channel span ratio estimation: %.3f\n', freqSpan/channSpan);

end



function fig10() % predictions on sweep variance and up/down asymmetry

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('pitchPredictionsSw.mat');

    swAsymm    = mean(asymm');
    swAsymmErr = std(asymm') ./ sqrt(3);

    f0 = [900, 1200, 1500];
    delta = linspace(-600, 600, 10);
    sSpace = 1:size(sp.f900d600, 2);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            condName = sprintf('f%.0fd%.0f', f0(i), delta(j));
            c{i}(j, :) = mean(sp.(condName));
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(j) = sum(sSpace .* fDist);
            sweepSigm(i, j) = (sum(sSpace.^2 .* fDist) - sweepPeak(j)^2);
        end
    end

    swFit = 1:10;

    load('deltaPitchDB.mat');

    f0  = avgFreqs{2};
    dSw = deltas{2};
    sw  = results;

    for i = 1:length(f0)
        FSw{i} = [];
        fSw{i} = [];
        for j = 1:length(subLab)
            FSw{i}  = [FSw{i};  sw{j}.F{i}(:, :) - f0(i)];
            fSw{i}  = [fSw{i};  sw{j}.F{i}(:, swFit) - f0(i)];
        end

        x = repmat(dSw(swFit), [length(subLab) * 2, 1]);
        X = repmat(dSw(swFit), [length(subLab) * 4, 1]);

    end

    fig = figure;
    subplot(121); hold off;
    for i = 1:length(subLab)
        for j = 1:length(f0) 
            swVar(length(f0) * (i - 1) + j, :) = std(sw{i}.F{j} - f0(j));
        end
    end
    yyaxis('left')
    x = {'-600', '', '', '-200', '', '', '200', '', '', '600'};
    distributionPlot(swVar, 'showMM', 4, 'xNames', x, ...
                     'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    ylabel('$std(\Delta p)$ (Hz)',  'Interpreter', 'LaTeX')
    xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
    ylim([0, 260])
    yticks(0:100:200);

    yyaxis('right')
    errorbar(1:10, mean(sweepSigm), std(sweepSigm) / sqrt(length(f0)), 'k', 'LineWidth', 1.5);
    ylabel('Var[channel]',  'Interpreter', 'LaTeX')
    ylim([.7, 6])
    yticks(1:1:5);

    swVarAbs = [swVar(:, 1:5); swVar(:, 10:-1:6)];
    swDeltaMat = repmat(linspace(600, 0, 5), [size(swVarAbs, 1), 1]);
    [swVarCorr, swVarP] = corr(swDeltaMat(:), swVarAbs(:), 'type', 'spearman');

    fprintf('Variance correlation (sweeps) = %.3f (p = %.2g))\n', swVarCorr, swVarP)

    subplot(122); hold off;
    fSwAll = [fSw{1}; fSw{2}; fSw{3}];
    asymmTmp = fSwAll(:, 1:5) + fSwAll(:, 10:-1:6);
    udSw = 0.25 * (asymmTmp(1:4:end,:) + asymmTmp(2:4:end,:) + ...
                   asymmTmp(3:4:end,:) + asymmTmp(4:4:end,:));
    yyaxis('left')
    x = {'600', '', '333', '', '66'};
    distributionPlot(udSw, 'showMM', 4, 'xNames', x, ...
                    'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (Hz)',  'Interpreter', 'LaTeX')
    yticks(-200:100:200);
    ylim([-200, 260])
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    uistack(zeroLine, 'bottom'); 
    yyaxis('right')
    errorbar(1:5, swAsymm, swAsymmErr, 'k', 'LineWidth', 1.5)
    ylim([-200, 260] / 10)
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (channels)',  'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2*plosCBwidths(2) 1.5];
    print(fig, 'fig2.svg', '-dsvg');

    legend({'model estimations', 'experimental data'}, 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2 * plosCBwidths(2) 1.75];
    print(fig, 'fig10.svg', '-dsvg');

end



function fig11() % parameter space exploration
    
    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('parameterSpace.mat'); % rr

    pars1 = {'$\tau^{memb}_e$', '$\tau^{memb}_e$', '$\Delta_{\omega s}$'};
    pars2 = {'$J_{NMDA}$', '$J_{NMDA}$', '$w_{\omega s}$ '};
    titles = {'adaptive $\tau_{pop}$', 'constant $\tau_{pop}$', ''};
    finalPars = {[0.05, 20], [], [0.05, 0.05]};

    fig = figure();
    for i = 1:length(r2)
        subplot(1,3,i)
        imagesc(r2{i}.par1, r2{i}.par2, r2{i}.R2', [0, 1])
        title(titles{i}, 'Interpreter', 'LaTeX')
        xlabel(pars1{i}, 'Interpreter', 'LaTeX')
        ylabel(pars2{i}, 'Interpreter', 'LaTeX')
        hold on
        if ~isempty(finalPars{i})
            plot(finalPars{i}(1), finalPars{i}(2), 'k+')
        end
    end    

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2*plosCBwidths(2) plosCBwidths(2)/3];
    print(fig, 'fig11.svg', '-dsvg');

    tpIx = find(r2{1}.par1 == finalPars{1}(1), 1);
    JnIx = find(r2{1}.par2 == finalPars{1}(2), 1);
    dfIx = find(r2{3}.par1 == finalPars{3}(1), 1);
    wfIx = find(r2{3}.par2 == finalPars{3}(2), 1);
    
    for i = 1:length(r2{3}.par1)
        for j = 1:length(r2{3}.par2)
            ix3(i, j) = (sqrt((i-dfIx)^2 + (j-wfIx)^2) < 5);
        end
    end                

    dwR2    = r2{3}.R2(ix3);
    jnR2    = r2{1}.R2(tpIx, (JnIx-5):(JnIx+5));
    finalR2 = r2{3}.R2(dfIx, wfIx);
    
    poolR2 = [jnR2(:); dwR2(:)];
    avgR2 = mean(poolR2);
    errR2 = std(poolR2) / sqrt(length(poolR2));

    meanR2  = sprintf('%.2f %s %.2f', avgR2, char(177), errR2);
    ratioR2 = sprintf('%.2f %s %.2f', avgR2/finalR2, char(177), errR2/finalR2);

    fprintf('Mean R2 across parametrization neighbourhood: %s\n', meanR2)
    fprintf('R2 ratio across parameter neighbourhood: %s\n', ratioR2)

end



function fig13() % Brady's stimuli

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('bradys.mat');
    sSpace = 1:size(bradyIImod{1,1}, 2);
    titles = {{'experiment II (down)',  'experiment II (up)'}, ...
              {'experiment III (down)', 'experiment III (up)'}};
    onsetsII  = 1000 * onsetsII;
    onsetsIII = 1000 * onsetsIII;

    fig = figure;
    for i = 1:2
        for j = 1:length(onsetsII)
            c{i}(j, :) = mean(bradyIImod{i, j});
            d{i}(j, :) = mean(bradyIIExp{i, j});
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            bradyIIPeak(j) = sum(sSpace .* fDist);
            bradyIISigma(j) = (sum(sSpace.^2 .* fDist) - bradyIIPeak(j)^2);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expPeakII(j) = sum(sSpace .* fDist);
        end

        subplot(1, 4, i)
        hold off
        imagesc(onsetsII, sSpace, c{i}');
        hold on
        errorbar(onsetsII, bradyIIPeak, bradyIISigma, 'LineWidth', 2, 'Color', 'black');

        plot(onsetsII, expPeakII, 's', 'Color', 'black', 'MarkerSize', 10);
        xticks(onsetsII)
        xlabel('transient onset (ms)', 'Interpreter', 'LaTeX')
        ylabel('cochlear channel index', 'Interpreter', 'LaTeX')
        ylim([30, 70])
        title(titles{1}{i}, 'Interpreter', 'LaTeX')
    end

    for i = 1:2
        for j = 1:length(onsetsIII)
            c{i}(j, :) = mean(bradyIIImod{i, j});
            d{i}(j, :) = mean(bradyIIIExp{i, j});
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            bradyIIIPeak(j) = sum(sSpace .* fDist);
            bradyIIISigma(j) = (sum(sSpace.^2 .* fDist) - bradyIIIPeak(j)^2);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expPeakIII(j) = sum(sSpace .* fDist);
        end

        subplot(1, 4, i+2)
        hold off
        imagesc(onsetsIII, sSpace, c{i}');
        hold on
        errorbar(onsetsIII, bradyIIIPeak, bradyIIISigma, 'LineWidth', 2, 'Color', 'black');
        plot(onsetsIII, expPeakIII, 's', 'Color', 'black', 'MarkerSize', 10);
        xticks(onsetsIII)
        xlabel('transient onset (ms)', 'Interpreter', 'LaTeX')
        ylabel('cochlear channel index', 'Interpreter', 'LaTeX')
        ylim([30, 70])
        title(titles{2}{i}, 'Interpreter', 'LaTeX')
    end

    legend('FM-feedback predictions', 'experimental results')


    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.5];
    print(fig, 'fig13.svg', '-dsvg');

    mo = [bradyIIPeak, bradyIIIPeak];
    ex = [expPeakII,   expPeakIII];
    R2 = 1 - sum((ex-mo).^2) / sum((ex - mean(ex)).^2);
    fprintf('FM-feedback model (brady): R2 = %.3f\n', R2);

end



function fig14() % train sweep pitch shift  

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    swFit = 1:6;

    load('deltaPitchSSDB.mat');

    f0  = avgFreqs{2};
    dSw = deltas{2};
    sw  = results;

    for i = 1:length(f0)
        FSw{i} = [];
        fSw{i} = [];
        for j = 1:length(subLab)
            FSw{i}  = [FSw{i};  sw{j}.F{i}(:, :) - f0(i)];
            fSw{i}  = [fSw{i};  sw{j}.F{i}(:, swFit) - f0(i)];
        end

        x = repmat(dSw(swFit), [length(subLab) * 2, 1]);
        X = repmat(dSw(swFit), [length(subLab) * 4, 1]);

        [p, s] = polyfit(X, fSw{i}, 1);
        fitErr = sqrt(diag(inv(s.R) * inv(s.R')) .* s.normr .^ 2 ./ s.df);
        fitResid = sum((fSw{i}(:) - polyval(p, X(:))).^2); 
        r2Sw(i)  = 1 - fitResid / ((length(fSw{i}(:)) - 1) * var(fSw{i}(:))); 
        fitSwm(i) = p(1);
        errSwm(i) = fitErr(1);
        fitSwn(i) = p(2);
        errSwn(i) = fitErr(2);

    end

    fig = figure;
    for i = 1:length(f0)
        subplot(1, length(f0), i); hold off;
        x = {'-333', '-200', '-66', '66', '200', '333'};
        distributionPlot(FSw{i}, 'showMM', 4, 'xNames', x, ...
                         'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
        title(['$\bar{f} = ' int2str(f0(i)) '$\,Hz'], 'Interpreter', 'LaTeX');
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('$\Delta p$ (Hz)', 'Interpreter', 'LaTeX')
        ylim([-333, 333])
        yticks(-300:100:300);
        zeroLine = refline(0, 0);
        zeroLine.Color = 'k';
        %fitLine = refline(fitSwm(i),fitSwn(i));
        m = fitSwm(i) * (dSw(3)- dSw(1)) / (3 - 1);
        n = fitSwm(i) * dSw(1) - m;
        fitLine = refline(m, n);
        fitLine.Color = .2 * [1 1 1];
        fitLine.LineWidth = 1;
        fitLine.LineStyle = '--';
        uistack(zeroLine, 'bottom'); 
    end

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig14.svg', '-dsvg');

end



function fig15() % heatmaps predictions train sweep pitch shift

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('pitchPredictionsTr.mat');
    
    f0 = [900, 1200, 1500];
    delta = linspace(-1000/3, 1000/3, 6);
    sSpace = 1:size(sp.f900d67, 2);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            condName = sprintf('f%.0fd%.0f', f0(i), delta(j));
            c{i}(j, :) = mean(sp.(condName));
            d{i}(j, :) = mean(ep.(condName));
        end
    end
    fig = figure;
    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(j) = sum(sSpace .* fDist);
            sweepSigm(j) = (sum(sSpace.^2 .* fDist) - sweepPeak(j)^2);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(j) = sum(sSpace .* fDist);
        end
        subplot(1, 3, i)
        hold off
        imagesc(delta, sSpace, c{i}');
        hold on
        errorbar(delta, sweepPeak, sweepSigm, 'LineWidth', 2, 'Color', 'black');

        plot(delta, expChannel, 's', 'Color', 'black', 'MarkerSize', 10);
        xticks([-600, -333.33, -66.66, 66.66, 333.33, 600])
        xticklabels([-600, -333, -66, 66, 333, 600])
        xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
        ylabel('cochlear channel index', 'Interpreter', 'LaTeX')
        ylim([25, 65])
        title(['$\bar{f} = ', sprintf('%d', f0(i)), '\,$Hz'], 'Interpreter', 'LaTeX')
    end

    legend('FM-feedback predictions', 'experimental results')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(3) 1.85];
    print(fig, 'fig15.svg', '-dsvg');

    for i = 1:length(f0)
        for j = 1:length(delta)
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(i,j) = sum(sSpace .* fDist);
            fDist = (exp(d{i}(j, :))) / sum(exp(d{i}(j, :)));
            expChannel(i,j) = sum(sSpace .* fDist);
        end
    end

    R2 = 1-sum((expChannel-sweepPeak).^2)/sum((expChannel-mean(expChannel)).^2);
    err2 = sqrt(mean((expChannel(:)-sweepPeak(:)).^2));

    fprintf('FM-feedback model (trains): R2 = %.3f, SME = %.3f\n', R2, err2);

end



function fig16() % model predictions up/down train asymmetry 

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('pitchPredictionsTr.mat');
    
    swAsymm    = mean(asymm');
    swAsymmErr = std(asymm') ./ sqrt(3);

    f0 = [900, 1200, 1500];
    delta = linspace(-1000/3, 1000/3, 6);
    sSpace = 1:size(sp.f900d67, 2);

    for i = 1:length(f0)
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
        for j = 1:length(delta)
            condName = sprintf('f%.0fd%.0f', f0(i), delta(j));
            c{i}(j, :) = mean(sp.(condName));
            fDist = (exp(c{i}(j, :))) / sum(exp(c{i}(j, :)));
            sweepPeak(j) = sum(sSpace .* fDist);
            sweepSigm(i, j) = (sum(sSpace.^2 .* fDist) - sweepPeak(j)^2);
        end
    end

    swFit = 1:10;
    swFit = 1:6;

    load('deltaPitchSSDB.mat');

    f0  = avgFreqs{1};
    dSw = deltas{1};

    for sIx = 1:length(subLab)
        [~, ix] = max(strcmp(subjects, subjects{sIx}));
        sw{sIx} = results{ix};
    end 

    for i = 1:length(f0)
        FSw{i} = [];
        fSw{i} = [];
        for j = 1:length(subLab)
            FSw{i}  = [FSw{i};  sw{j}.F{i}(:, :) - f0(i)];
            fSw{i}  = [fSw{i};  sw{j}.F{i}(:, swFit) - f0(i)];
        end

        x = repmat(dSw(swFit), [length(subLab) * 2, 1]);
        X = repmat(dSw(swFit), [length(subLab) * 4, 1]);

    end

    fig = figure;
    subplot(121); hold off;
    for i = 1:length(subLab)
        for j = 1:length(f0) 
            swVar(length(f0) * (i - 1) + j, :) = std(sw{i}.F{j} - f0(j));
        end
    end
    yyaxis('left')
    x = {'-333', '-200', '-66', '66', '200', '333'};
    distributionPlot(swVar, 'showMM', 4, 'xNames', x, ...
                     'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    ylabel('$std(\Delta p)$ (Hz)',  'Interpreter', 'LaTeX')
    xlabel('$\Delta f$ (Hz)', 'Interpreter', 'LaTeX')
    ylim([0, 220])
    yticks(0:100:200);

    yyaxis('right')
    errorbar(1:6, mean(sweepSigm), std(sweepSigm) / sqrt(length(f0)), 'k', 'LineWidth', 1.5);
    ylabel('Var[channel]',  'Interpreter', 'LaTeX')
    ylim([1.1, 12])
    yticks(2:2:12);


    swVarAbs = [swVar(:, 1:3); swVar(:, 6:-1:4)];
    swDeltaMat = repmat(linspace(333, 0, 3), [size(swVarAbs, 1), 1]);
    [swVarCorr, swVarP] = corr(swDeltaMat(:), swVarAbs(:), 'type', 'spearman');

    fprintf('Variance correlation (trains) = %.3f (p = %.2g))\n', swVarCorr, swVarP)

    subplot(122); hold off;
    fSwAll = [fSw{1}; fSw{2}; fSw{3}];
    asymmTmp = fSwAll(:, 1:3) + fSwAll(:, 6:-1:4);
    udSw = 0.25 * (asymmTmp(1:4:end,:) + asymmTmp(2:4:end,:) + ...
                   asymmTmp(3:4:end,:) + asymmTmp(4:4:end,:));
    yyaxis('left')  
    x = {'333', '200', '66'};
    distributionPlot(udSw, 'showMM', 4, 'xNames', x, ...
                    'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (Hz)',  'Interpreter', 'LaTeX')
    ylim([-290, 250])
    yticks(-200:100:200);
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    uistack(zeroLine, 'bottom'); 
    
    yyaxis('right')
    errorbar(1:3, swAsymm, swAsymmErr, 'k', 'LineWidth', 1.5)
    ylim([-290, 250]/10)
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (channels)',  'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2*plosCBwidths(2) 1.5];
    print(fig, 'fig2.svg', '-dsvg');

    legend({'model estimations', 'experimental data'}, 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2 * plosCBwidths(2) 1.75];
    print(fig, 'fig16.svg', '-dsvg');

end



function fig17() % Stimulus example 

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    fs        = 48000;
    durations = [0.005, 0.04, 0.005]; 
    fAvg      = 1200;
    Delta     = 200;
    T         = 0.005;

    tSw = 1000 * linspace(0, sum(durations), fs * sum(durations));
    tTr = 1000 * linspace(0, sum(durations) + 4 * durations(2), ...
                          fs * (sum(durations) + 4 * durations(2)));
    N  = durations * fs;
    f0 = fAvg + 0.5 * Delta;
    f1 = fAvg - 0.5 * Delta;

    f      = 1 ./ linspace(1 / f0, 1 / f1, N(2));
    fSweep = [f0 * ones([1, N(1)]), f, f1 * ones([1, N(3)])];
    phase  = 2 * pi * cumsum(fSweep) / fs;
    sSweep = 1 / length(f0) * sin(phase);

    fTrain = [f0 * ones([1, N(1)]), repmat(f, [1, 5]), f1 * ones([1, N(3)])];
    phase  = 2 * pi * cumsum(fTrain) / fs;
    sTrain = 1 / length(f0) * sin(phase);

    hamm = hamming(2 * T * fs);
    modd = ones(size(sSweep));
    modd(1:(T * fs)) = hamm(1:(T * fs));
    modd((length(sSweep) - (T*fs) + 1):end) = hamm(((T*fs)+1):end);
    sSweep = sSweep .* modd;
    modd = ones(size(sTrain));
    modd(1:(T * fs)) = hamm(1:(T * fs));
    modd((length(sTrain) - (T*fs) + 1):end) = hamm(((T*fs)+1):end);
    sTrain = sTrain .* modd;

    fig = figure;
    subplot(4,6,1); 
    plot(tSw, sSweep, 'k', 'LineWidth', 0.5); 
    title('sweep', 'Interpreter', 'LaTeX')
    ylabel('amplitude', 'Interpreter', 'LaTeX'); 
    xticks([5, 45])
    xticklabels({})
    ylim([-1.1, 1.1])
    xlim([0, 1000 * (sum(durations))])
    
    subplot(4,6,2:6); 
    plot(tTr, sTrain, 'k', 'LineWidth', 0.5);
    title('sweep train', 'Interpreter', 'LaTeX')
    xticks([5, 45, 85, 125, 165, 205])
    xticklabels({})
    yticklabels({})
    ylim([-1.1, 1.1])
    xlim([0, 1000 * (sum(durations) + 4 * durations(2))])


    subplot(4,6,[7, 13, 19]); 
    plot(tSw, fSweep / 1000, 'k', 'LineWidth', 2); 
    ylabel('frequency (kHz)', 'Interpreter', 'LaTeX')
    ylim([1.090, 1.310])
    xlim([0, 1000 * (sum(durations))])
    xticks([5, 45])
    xlabel('time (ms)', 'Interpreter', 'LaTeX')

    subplot(4,6,[8:12, 14:18, 20:24]); 
    plot(tTr, fTrain / 1000, 'k', 'LineWidth', 2); 
    ylim([1.090, 1.310])
    xlim([0, 1000 * (sum(durations) + 4 * durations(2))])
    xticks([5, 45, 85, 125, 165, 205])
    xlabel('time (ms)', 'Interpreter', 'LaTeX')
    yticklabels({})


    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.25*plosCBwidths(3) 1.85];
    print(fig, 'fig17.svg', '-dsvg');

end



function fig18() % Connectivity matrices

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    load('./connectivity.mat');

    fig = figure;
    interval = 1:25;

    subplot(2, 3, 1)
    imagesc(interval, interval, downei(interval, interval))
    title('$\omega^{ei}$', 'Interpreter', 'LaTeX')
    
    subplot(2, 3, 4)
    imagesc(interval, interval, downie(interval, interval))
    title('$\omega^{ie}$', 'Interpreter', 'LaTeX')
    ylabel('target population', 'Interpreter', 'LaTeX')
    
    subplot(2, 3, 2)
    imagesc(interval, interval, fup(interval, interval))
    title('$\omega^{f\uparrow}$', 'Interpreter', 'LaTeX')
    
    subplot(2, 3, 5)
    imagesc(interval, interval, fdown(interval, interval))
    title('$\omega^{f\downarrow}$', 'Interpreter', 'LaTeX')
    xlabel('source population', 'Interpreter', 'LaTeX')

    subplot(2, 3, 3)
    imagesc(interval, interval, upf(interval, interval))
    title('$\omega^{\uparrow f}$', 'Interpreter', 'LaTeX')

    subplot(2, 3, 6)
    imagesc(interval, interval, downf(interval, interval))
    title('$\omega^{\downarrow f}$', 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4.5 3];
    print(fig, 'fig18.svg', '-dsvg');

    fig = figure;
    subplot(131)
    imagesc([0, 1; 1, 0])
    c = colorbar();
    c.Label.String = 'activity distribution';
    c.Label.Interpreter = 'LaTeX';
    c.Ticks = [0, 1];

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 4.5 1.7];
    print(fig, 'fig18colorbar.svg', '-dsvg');

end



function figS1() % Effect of presentation order on Delta p

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    swFit = 1:10;
    trFit = 1:6;

    matDBsw = load('deltaPitchDB.mat');
    subLab = matDBsw.subLab;
    f0  = matDBsw.avgFreqs{2};
    dSw = matDBsw.deltas{2};
    sw = matDBsw.results;

    for i = 1:length(f0)
        f1Sw{i} = []; f2Sw{i} = [];
        for j = 1:length(subLab)
            f1Sw{i} = [f1Sw{i}; sw{j}.F{i}(1:2, swFit) - f0(i)];
            f2Sw{i} = [f2Sw{i}; sw{j}.F{i}(3:4, swFit) - f0(i)];
        end
    end
    f1SwAll = [f1Sw{1}; f1Sw{2}; f1Sw{3}];
    f2SwAll = [f2Sw{1}; f2Sw{2}; f2Sw{3}];
    asymmTmp = f1SwAll - f2SwAll;
    lrSw = 0.5 * (asymmTmp(1:2:end, :) + asymmTmp(2:2:end, :));
    
    fig = figure;
    subplot(121); hold off;
    x = {'-600', '', '-333', '', '-66', '66', '', '333', '', '600'};
    distributionPlot(lrSw, 'showMM', 4, 'xNames', x, ...
                    'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (Hz)',  'Interpreter', 'LaTeX')
    yticks(-200:100:200);
    ylim([-200, 260])
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    uistack(zeroLine, 'bottom'); 
    title('sweeps', 'Interpreter', 'LaTeX')

    matDBsw = load('deltaPitchSSDB.mat');

    f0  = matDBsw.avgFreqs{2};
    dSw = matDBsw.deltas{2};
    sw  = matDBsw.results;

    for i = 1:length(f0)
        f1Sw{i} = []; f2Sw{i} = [];
        for j = 1:length(subLab)
            f1Sw{i} = [f1Sw{i}; sw{j}.F{i}(1:2, trFit) - f0(i)];
            f2Sw{i} = [f2Sw{i}; sw{j}.F{i}(3:4, trFit) - f0(i)];
        end
    end
    f1SwAll = [f1Sw{1}; f1Sw{2}; f1Sw{3}];
    f2SwAll = [f2Sw{1}; f2Sw{2}; f2Sw{3}];
    asymmTmp = f1SwAll - f2SwAll;
    lrSw = 0.5 * (asymmTmp(1:2:end, :) + asymmTmp(2:2:end, :));
    
    subplot(122); hold off;
    x = {'-333', '-200', '-66', '66', '200', '300'};
    distributionPlot(lrSw, 'showMM', 4, 'xNames', x, ...
                    'histOpt', 1, 'color', plosCBcolors(1, :), 'distWidth', 1)
    xlabel('$|\Delta f|$ (Hz)', 'Interpreter', 'LaTeX')
    ylabel('asymm$^{\uparrow\downarrow}_{|\Delta f|}$ (Hz)',  'Interpreter', 'LaTeX')
    ylim([-290, 250])
    yticks(-200:100:200);
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    uistack(zeroLine, 'bottom'); 
    title('sweep trains', 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.35*plosCBwidths(2) 1.5];
    print(fig, 'figS1.svg', '-dsvg');


end



function figS2() % Slope of the linear fits per subject

    [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling();

    matDBsw = load('deltaPitchDB.mat');
    matDBtr = load('deltaPitchSSDB.mat');

    subjects = matDBsw.subjects;
    subLab   = matDBsw.subLab;

    f0  = matDBsw.avgFreqs{2};
    dSw = matDBsw.deltas{2};
    dTr = matDBtr.deltas{2};

    sw = matDBsw.results;
    tr = matDBtr.results;

    % Slopes per subject and centre frequency
    for j = 1:length(subLab);
        for i = 1:length(f0)
            [p, s] = polyfit(repmat(dSw, [4, 1]), sw{j}.F{i} - f0(i), 1);
            fitErr = sqrt(diag(inv(s.R) * inv(s.R')) .* s.normr .^ 2 ./ s.df);
            mSwErr(i, j)  = fitErr(1);
            mSw(i, j) = p(1);
            [p, s] = polyfit(repmat(dTr, [4, 1]), tr{j}.F{i} - f0(i), 1);
            fitErr = sqrt(diag(inv(s.R) * inv(s.R')) .* s.normr .^ 2 ./ s.df);
            mTrErr(i, j)  = fitErr(1);
            mTr(i, j) = p(1);
        end
    end 

    fig = figure;
    subplot(121); hold off
    styles = ['so^'];
    for i = 1:length(f0)
        errorbar(1:length(subLab), mSw(i, :), mSwErr(i, :), 'Color', plosCBcolors(i, :), 'Marker', styles(i));
        hold on
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
    end
    xlim([0.5, 8.5])
    xlabel('subjects', 'Interpreter', 'LaTeX')
    ylabel('slope of the linear fit between $f_{perceived}$ and $\Delta f$', 'Interpreter', 'LaTeX')
    ylim([-0.35, 0.65])
    xticks(1:8)
    xticklabels(subLab)
    l = legend(freqLegend, 'Interpreter', 'LaTeX', 'Location', 'SouthEast');
    l.AutoUpdate = 'off';
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    zeroLine.LineWidth = 0.8;
    uistack(zeroLine, 'bottom'); 
    title('single sweeps', 'Interpreter', 'LaTeX')

    subplot(122); hold off
    styles = ['so^'];
    for i = 1:length(f0)
        errorbar(1:length(subLab), mTr(i, :), mTrErr(i, :), 'Color', plosCBcolors(i, :), 'Marker', styles(i));
        hold on
        freqLegend{i} = sprintf('$\\bar{f} = %.0f$ Hz', f0(i));
    end
    xlim([0.5, 8.5])
    xlabel('subjects', 'Interpreter', 'LaTeX')
    ylabel('slope of the linear fit between $f_{perceived}$ and $\Delta f$', 'Interpreter', 'LaTeX')
    ylim([-0.35, 0.65])
    xticks(1:8)
    xticklabels(subLab)
    zeroLine = refline(0, 0);
    zeroLine.Color = 'k';
    zeroLine.LineWidth = 0.8;
    uistack(zeroLine, 'bottom');
    title('sweep trains', 'Interpreter', 'LaTeX')

    set(findall(gcf,'-property','FontSize'),'FontSize',plosCBfontsize)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 1.2*plosCBwidths(2) 1.5];
    print(fig, 'figS2.svg', '-dsvg');

end



% Auxiliary functions

function figureCanvasForPlos()

    plosCBwidths = styling();

    figTemplate = figure;
    for i = 1:length(plosCBwidths)
        plot([0, 1], [0, 1])
        figTemplate.PaperUnits = 'inches';
        figTemplate.PaperPosition = [0 0 plosCBwidths(i) 10];
        print(figTemplate, sprintf('figTemplate%d.svg', i), '-dsvg');
    end

end

    
function [plosCBwidths, plosCBfontsize, plosCBcolors, styles] = styling()

    plosCBwidths   = [2.63, 5.2, 7.5]; % inches
    plosCBfontsize = 8; 
    plosCBcolors   = [0, 142, 255; 255, 69, 0; 13, 185, 0] / 255;
    styles = ['so^'];

end
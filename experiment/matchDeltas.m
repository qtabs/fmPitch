function matchDeltas(resume)

    if nargin < 1, resume = 0; end

    
    pars  = getParameters();
    mess  = getMessages(pars);

    if resume == 1
        cache   = load('cache.mat');
        blcks   = cache.blcks;
        subject = cache.subject;
        r       = cache.r;
        n0      = cache.n + 1;
        clear cache;
        fprintf('%s: %s; %s: %d\n', mess.subj, subject, mess.lastBlock, n0-1);
        x = input(sprintf('%s  ', mess.repeatTr), 's');
        if strncmpi(x, 'y', 1) || strncmpi(x, 'j', 1) || strncmp(x, '1', 1)
            training = 1;
        else
            training = 0;
        end
        input(sprintf(' -- %s...', mess.pressEnter), 's');
    else
        subject = input(sprintf('%s: ', mess.subj), 's');
        blcks   = generateBlocks(pars);
        n0      = 1;
        training = 1;
    end

    if not(strcmp(subject, 'lazy')) && (training == 1)
        fprintf('\n\n---------------------------------------- ');
        fprintf(mess.training);
        fprintf(' ----------------------------------------\n');
        trainSucc = playTraining(pars, mess);
        if trainSucc
            fprintf('\n---------------------------------------- ');
            fprintf (mess.endTrain);
            fprintf(' ----------------------------------------\n');
        else
            return
        end
            
    end

    if n0 == 1   
          
        fprintf('\n\n---------------------------------------- ');
        fprintf('%s %s', mess.test, mess.block);
        fprintf(' ----------------------------------------\n');
        
        [r{1}, testSucc] = playTest(blcks{1}, pars, mess);

        if testSucc
            n = 1; n0 = 2;
            save('cache.mat', 'r', 'blcks', 'subject', 'n');
        
            fprintf('\n---------------------------------------- ');
            fprintf (mess.endTest);
            fprintf(' ----------------------------------------\n');
            
            fprintf('\n%s ', mess.break);
            fprintf(mess.continue)
            pause();       
        else
            return
        end
        
    end

    for n = n0:length(blcks)
        
        fprintf('\n\n---------------------------------------- ');
        fprintf('%s %d/%d', mess.block, n, length(blcks));
        fprintf(' ----------------------------------------\n');
        
        r{n} = playBlock(blcks{n}, pars, mess);

        save('cache.mat', 'r', 'blcks', 'subject', 'n');
        
        fprintf('---------------------------------------- ');
        fprintf(mess.endBlock);
        fprintf(' ----------------------------------------\n');

        if n < length(blcks)
            fprintf('\n%s ', mess.break);
            fprintf(mess.continue)
            pause();
        end
    end

    saveResults(r, subject, pars);
    fprintf('%s\n\n', mess.done)

end


function blcks = generateBlocks(pars)

    [deltas, avgFreqs, testCombs] = getSpaces(pars);

    blck = combvec(deltas, avgFreqs)';

    testBlk = [];
    for i = 1:size(testCombs, 1)
        for j = 1:size(blck, 1)
            if sum(abs(blck(j, :) - testCombs(i, :))) < 1
                testBlk(i, :) = blck(j, :);
                blck = [blck(1:(j-1), :); blck((j+1):end, :)];
                break
            end
        end
    end

    dirs = [-1 * ones([floor(size(blck, 1) / 2), 1]); ...
                 ones([ ceil(size(blck, 1) / 2), 1])];
    dirs = dirs(randperm(length(dirs)));

    for i = 1:(2 * pars.N)
        dirs = -1 * dirs;
        blck(:, 3) = dirs;
        blcks{i + 1} = blck(randperm(size(blck, 1)), :);
    end

    dirs = [-1 * ones([floor(size(testBlk, 1) / 2), 1]); ...
                 ones([ ceil(size(testBlk, 1) / 2), 1])];
    dirs = dirs(randperm(length(dirs)));
    for i = 1:(2 * pars.N)
        dirs = -1 * dirs;
        testBlk(:, 3) = dirs;
        sl = (1 + (i - 1) * size(testBlk, 1)):(i * size(testBlk, 1));
        blcks{1}(sl, :) = testBlk(randperm(size(testBlk, 1)), :);
    end

end



function r = playBlock(blck, pars, mess)

    for i = 1:size(blck, 1)
        fprintf(' -- %s %d/%d --\n', mess.sound, i, size(blck, 1)); 
        [F(i), T(i), P{i}] = graduateRef(blck(i, :), pars, mess);
        fprintf('\n');
        save('blockCache.mat', 'F', 'T', 'P', 'blck', 'i');
        pause(0.5);
    end

    blck(:, 4) = 1:size(blck, 1);
    r.F        = F;
    r.T        = T;
    r.P        = P;
    r.specs    = blck;

end



function trainSucc = playTraining(pars, mess)

    thresh = 1.1 * pars.epsilon;
    err    = thresh + 1; 
    trials = 0;
    [~, avgFreqs] = getSpaces(pars);

    if strcmp(pars.fSpace, 'inverse')
        blck = repmat(combvec(1, avgFreqs)', [2, 1]);
    else
        blck = repmat(combvec(0, avgFreqs)', [2, 1]);
    end

    blck = blck(randperm(size(blck, 1)), :);
    dirs = [-1 * ones([floor(size(blck, 1) / 2), 1]); ...
                 ones([ ceil(size(blck, 1) / 2), 1])];
    dirs = dirs(randperm(length(dirs)));

    blck(:, 2) = blck(:, 2) .* (1 + 0.2 * (randn(size(blck(:, 2)))));
    blck(:, 2) = 50 * round(blck(:, 2) / 50);
    blck(:, 3) = dirs;

    for i = 1:size(blck, 1)
        fprintf(' -- %s %d/%d --\n', mess.sound, i, size(blck, 1)); 
        [F(i)] = graduateRef(blck(i, :), pars, mess, 1);
        fprintf('\n');
        pause(0.5);
    end

    while (err > thresh) && (trials < 6)
        fprintf('%s...\n', mess.trContinue);

        blck = blck(randperm(size(blck, 1)), :);
        dirs = [-1 * ones([floor(size(blck, 1) / 2), 1]); ...
                     ones([ ceil(size(blck, 1) / 2), 1])];
        dirs = dirs(randperm(length(dirs)));

        blck(:, 2) = blck(:, 2) .* (1 - 0.1 * abs(randn(size(blck(:, 2)))));
        blck(:, 2) = 50 * round(blck(:, 2) / 50);
        blck(:, 3) = dirs;

        for i = 1:size(blck, 1)
            fprintf(' -- %s %d/%d --\n', mess.sound, i, size(blck, 1)); 
            [F(i)] = graduateRef(blck(i, :), pars, mess, 0);
            fprintf('\n');
            pause(0.5);
        end

        trials = trials + 1;
        err    = sum(abs(F' - blck(:, 2)));
    end

    if err < thresh
        fprintf('%s %s...', mess.trainDone, mess.pressEnter);
        trainSucc = true;
        pause();
    else 
        fprintf(' -- %s --\n', mess.trFail);
        trainSucc = false;
    end

end



function [r, testSucc] = playTest(blck, pars, mess)

    [~, ~, testCombs] = getSpaces(pars);
    r = playBlock(blck, pars, mess);

    for i = 1:size(testCombs, 1)
        k = 0;
        for j = 1:size(blck, 1)
            if sum(abs(blck(j, 1:2) - testCombs(i, :))) < 1
                k = k + 1;
                F(k, i) = r.F(j);
            end 
        end
    end

    for i = 1:size(F, 1)
        A(:, i) = std(F);
    end
    constist = 1 - (2 * (sqrt(mean(A(:).^2))) / 175);
    fprintf('[Constistency = %.1f]\n', constist)

    if constist > 0.50
        fprintf('%s %s...', mess.trainDone, mess.pressEnter);
        testSucc = true;
        pause();
    else 
        fprintf(' -- %s --\n', mess.ttFail);
        testSucc = false;
    end

end



function [fRef, time, gradPath] = graduateRef(specs, pars, mess, hints)

    if nargin < 4, hints = -1; end;

    if strcmp(pars.fSpace, 'inverse')
        f0 = specs(2) * specs(1);
        f1 = specs(2) / specs(1);
    else
        f0 = specs(2) - 0.5 * specs(1);
        f1 = specs(2) + 0.5 * specs(1);
    end

    fRef = specs(2) * min(1.5, max(0.5, 1 + 0.15 * randn()));
    fRef = pars.epsilon * round(fRef / pars.epsilon);

    modif = -1;
    tTrack = tic; gradPath = [];

    while not(modif == 0)

        if specs(3) < 0
            playStim(f0, f1, pars);
            pause(pars.ISI)
            playStim(fRef, fRef, pars);
        else
            playStim(fRef, fRef, pars);
            pause(pars.ISI)
            playStim(f0, f1, pars);
        end

        if hints > 0
            if (fRef > f0) & (specs(3) < 0)
                hint = mess.higher;
            elseif (fRef > f0) & (specs(3) > 0)
                hint = mess.lower;
            elseif (fRef < f0) & (specs(3) < 0)
                hint = mess.lower;
            elseif (fRef < f0) & (specs(3) > 0)
                hint = mess.higher;
            else
                hint = mess.same;
            end
            hint = sprintf('(%s %s)', mess.hint, hint);
            modif = getFeedback(pars, mess, hint);
        else
            modif = getFeedback(pars, mess);
        end
        
        if modif == 2
            fRef = fRef + specs(3) * pars.epsilon;
        elseif modif == 1
            fRef = fRef - specs(3) * pars.epsilon;
        end

        gradPath = [gradPath, modif];

    end

    if hints >= 0
        if fRef == specs(2)
            fprintf('%s!\n', mess.correct)
        else
            fprintf('%s\n', mess.wrong)
        end
    end

    time = toc(tTrack);

end



function x = getFeedback(pars, mess, hint)

    if nargin < 3
        feedbackMess = mess.feedback;
    else
        feedbackMess = [mess.feedback, '\n', hint, ' '];
    end

    x  = input(feedbackMess, 's');

    % Checking answer
    [x, status] = str2num(x);

    if isempty(x)
        x = -1;
    elseif size(x, 1) > 1 || size(x, 2) > 1
        x = -2;
        fprintf('%s %s.\n', mess.errMess, mess.singNumb);
        pause(0.5);
    elseif x < 0 || x > 2
        x = -2;
        fprintf('%s %s.\n', mess.errMess, mess.integer);
        pause(0.5);
    end

end



function playStim(f0, f1, pars)
   
    N = pars.durations * pars.fs;
    stim = zeros([1, sum(N)]);
    
    for i = 1:length(f0)
        if strcmp(pars.stimSpace, 'inverse')
            fsweep = 1 ./ linspace(1 / f0(i), 1 / f1(i), N(2));
        else
            fsweep = linspace(f0(i), f1(i), N(2)); 
        end
        f = [f0(i) * ones([1, N(1)]), fsweep, f1(i) * ones([1, N(3)])];
        phase = 2 * pi * cumsum(f) / pars.fs;
        stim = stim + (1 / length(f0)) * sin(phase);
    end

    stim = ramps(stim, 0.005, pars.fs);

    %%%%%%%%%%%%%%%
    % c = clock();
    % audiowrite(sprintf('sound%02.0f%2.0f.wav', c(5), c(6)), stim, pars.fs);
    %%%%%%%%%%%%%%%

    sound(stim, pars.fs);
    pause(sum(pars.durations) + 0.01);

end



function stim = ramps(stim, T, fs)

    if nargin < 3, fs = 48000; end;
    if nargin < 2, T = 0.005; end;

    N = T * fs;
    
    hamm = hamming(2 * N);
    modd = ones(size(stim));

    modd(1:N)                      = hamm(1:N);
    modd((length(stim) - N + 1):end) = hamm((N+1):end);

    stim = [zeros([1, 0.05*fs]), stim .* modd, zeros([1, 0.05*fs])];

end



function saveResults(r, subject, pars)

    [deltas, avgFreqs] = getSpaces(pars);

    for i = 1:length(avgFreqs)
        for j = 1:length(deltas)    
            m = 0;
            for direction = [-1, 1]
                for n = 1:length(r)
                    for k = 1:length(r{n}.F)
                        if (r{n}.specs(k, 2) == avgFreqs(i)) && ...
                           (abs(r{n}.specs(k, 1) - deltas(j)) < 0.1) && ...
                           (r{n}.specs(k, 3) == direction)
                            m = m + 1;
                            results.F{i}(m, j)  = r{n}.F(k);
                            results.T{i}(m, j)  = r{n}.T(k);
                            results.P{i}{m, j}  = r{n}.P{k};
                            results.ix{i}(m, j) = r{n}.specs(k, 4);
                            results.d{i}(m, j)  = r{n}.specs(k, 3);
                        end
                    end
                end
            end
        end
    end

    
    for i = 1:length(avgFreqs)
        A(:, i) = std(results.F{i}(:, 4:7));
    end

    constistency = 1 - (2 * (sqrt(mean(A(:).^2))) / 175);
    fprintf('[Constistency = %.1f]\n\n', constistency);

    if exist('deltaPitchDB.mat', 'file') == 0
        matDB.subjects = {};
        save('deltaPitchDB', '-struct', 'matDB');
    end

    matDB    = load('deltaPitchDB.mat');
    lenMatDB = length(matDB.subjects);

    matDB.subjects{lenMatDB + 1} = subject;
    matDB.results{lenMatDB + 1}  = results;
    matDB.avgFreqs{lenMatDB + 1} = avgFreqs;
    matDB.deltas{lenMatDB +1 }   = deltas;

    save('deltaPitchDB', '-struct', 'matDB');

end



function pars = getParameters()

    pars.fs        = 48000; % Herz
    pars.durations = [0.005, 0.04, 0.005]; % seconds
    pars.epsilon   = 25;    % Herz
    pars.N         = 2;
    pars.ISI       = 0.1; % ms
    pars.fSpace    = 'linear';
    pars.stimSpace = 'inverse';
    pars.language  = 'deutsch';

end



function mess = getMessages(pars)

    if strcmp(pars.language, 'english')

        mess.feedback = 'Is the last sound higher (2), lower (1), or of the ';
        mess.feedback = [mess.feedback, 'same pitch (0) as the first one? '];
        mess.feedback = [mess.feedback, 'Press ENTER to repeat. '];
        mess.lastBlock  = 'last block in cache';
        mess.repeatTr   = 'Do you want to repeat the training session?';
        mess.pressEnter = 'Press ENTER to continue';
        mess.subj       = 'Subject reference';
        mess.training   = 'TRAINING SESSION';
        mess.endTrain   = 'end of training';
        mess.block      = 'BLOCK';
        mess.endBlock   = 'end of block';
        mess.break      = 'Please, take a break!';
        mess.continue   = 'Press ENTER to continue when ready :)';
        mess.done       = 'Done! Thanks! <3 <3';
        mess.sound      = 'Sound';
        mess.trContinue = 'Continuing training';
        mess.trainDone  = 'Training completed!';
        mess.errMess    = 'Answer not recognised! Repeating sound...';
        mess.singNumb   = 'please, next time enter a single number';
        mess.integer    = 'please, next time enter 0, 1, or 2';
        mess.trFail     = 'Training unsuccessful';
        mess.test       = 'TEST';
        mess.ttFail     = 'Test unsuccessful';
        mess.endTest    = 'end of test';


        mess.hint       = 'hint: the second tone is now';
        mess.same       = 'as high/low as the first one';
        mess.higher     = 'higher';
        mess.lower      = 'lower';
        mess.correct    = 'Correct';
        mess.wrong      = 'Wrong';

    elseif strcmp(pars.language, 'deutsch')

        mess.feedback = 'Druecken: 1 (zweiter Ton tiefer), 2 (zweiter hoeher),';
        mess.feedback = [mess.feedback,' 0 (gleiche Tonhoehe), oder ENTER'];
        mess.feedback = [mess.feedback,' um zu wiederholen. '];
        mess.lastBlock  = 'letzter Block im Zwischenspeicher';
        mess.repeatTr   = 'Moechten Sie das Training wiederholen?';
        mess.pressEnter = 'Druecken Sie ENTER um fortzufahren';
        mess.subj       = 'Proband(in) ref';
        mess.training   = 'TRAININGSDURCHLAUF';
        mess.endTrain   = 'Ende des Trainings';
        mess.block      = 'BLOCK';
        mess.endBlock   = 'Ende des Blocks';
        mess.break      = 'Bitte, machen Sie eine Pause!';
        mess.continue   = 'Druecken Sie ENTER um fortzufahren, wenn Sie ';
        mess.continue   =  [mess.continue, 'bereit sind.'];
        mess.done       = 'Fertig! Vielen Dank!';
        mess.sound      = 'Ton';
        mess.trContinue = 'Fortfahren mit dem Training';
        mess.trainDone  = 'Training abgeschlossen!';
        mess.errMess    = 'Antwort nicht erkannt, Ton wird wiederholt...';
        mess.singNumb   = 'bitte druecken Sie beim naechsten Durchlauf eine ';
        mess.singNumb   = [mess.singNumb, 'einzelne Nummer'];
        mess.integer    = 'bitte druecken Sie beim naechsten Durchlauf 0, 1 ';
        mess.integer    = [mess.integer, 'oder 2.'];
        mess.trFail     = 'Training nicht erfolgreich';
        mess.test       = 'TEST';
        mess.ttFail     = 'Test nicht erfolgreich';
        mess.endTest    = 'Ende des Test';

        mess.hint       = 'Hinweis: der zweite Ton ist nun';
        mess.same       = 'so hoch/tief wie der erste';
        mess.higher     = 'hoeher';
        mess.lower      = 'tiefer';
        mess.correct    = 'Richtig';
        mess.wrong      = 'Falsch';

    end

end



function [deltas, avgFreqs, testCombs] = getSpaces(pars)
    
    % Hardcoded parameters (do not touch! The entire DB depends on it!)
    if strcmp(pars.fSpace, 'inverse')
        deltas = linspace(1/1.5, 1.5, 10);
    else
        deltas = linspace(-600, 600, 10);
    end

    avgFreqs  = [900, 1200, 1500];
    testCombs = [67, 900; -200, 1200; -67, 1500];

end



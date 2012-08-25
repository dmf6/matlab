function [frequencies, spikephases, spiketimes, isis, avgisis] = spike_phase_generator_vclamp(filename, pre)


    %out = load('2012_07_17_22_31_47.asc');
    %out = load('2012_07_30_20_07_18.asc');
    %out = load('2012_07_30_20_09_56.asc');
    %out = load('2012_07_30_20_42_24.asc');
    %out = load('2012_07_30_20_45_04.asc');
    %out = load('2012_07_30_21_04_18.asc');
    %out = load('2012_08_01_17_35_38.asc');
    %out = load('2012_08_01_19_05_57.asc');
    %out = load('2012_08_01_20_11_53.asc');
    %out = load('2012_08_01_20_27_14.asc');
    out = load(filename);

    t = out(:, 1)-pre;
    vv = out(:, 2);
    ii = out(:, 3);
  
    %chop off pre recording, start from 0
    tidx = find(t > -0.01);
    ti = tidx(1);
    
     t=t(ti:length(tidx));
     vv=vv(ti:length(tidx));
     ii=ii(ti:length(tidx));
%     
     N = length(t);
     dt = t(2) - t(1);
% 
     T = N*dt; % Number of samples
     fs = N/(T/1000); %Sample frequency divide by 1000 to convert to Hz
% 
     fc = 100; % cut-off frequency
% 
     [b, a] = butter(2, fc/fs, 'high');
     icut = filter(b, a, ii);
% 
     figure(1)
     subplot(3, 1, 1);
     plot(t, ii, '-r');
     axis([0 t(length(t)) -0.2 0.2])
% 
     subplot(3, 1, 2);
     plot(t, icut, '-b');
     hold on 
%     % define threshold for a spike
% 
     i_thr = mean(icut) + 2*std(icut);
     ithr = 0:50:t(length(t));
     ithr(:) = i_thr;
     t_new = 0:50:t(length(t));
     plot(t_new, ithr, '--k', 'Linewidth', 3);
     axis([0 t(length(t)) -0.2 0.2])
%     
% 
     subplot(3, 1, 3)
% 
     plot(t, vv, '-g');
     axis([0 t(length(t)) -0.7 -0.2])
% 
% % % 
      v_thr = -0.55; 
      a = find(vv > v_thr);
     cycle_times = t(a(1));
     last_cross = cycle_times;
     cyclecount=1;
     flag = 0;
% 
     for i=a(2):length(ii)
%         %noise is 30ms, that is we could hit the same current reading in a time
%         %window of 30ms. If that happens ignore
             if (vv(i) > v_thr && flag == 1)
                 if (abs(t(i)-last_cross) > 20)
                     cycle_times = [cycle_times t(i)];
                     cyclecount=cyclecount + 1;  
                 end
                 last_cross = t(i);
                 flag=0;
             elseif(vv(i) < v_thr && flag == 0)
                  last_cross = t(i);  
                  flag=1;
             end
     end
% % 
% %     %current trace is rather noisy. Consider using ME1 for current injection
% %     %(for more resolution) and recoring the injected current not the measured
% %     %current (inputs on back of amp not outputs)
% % 
     cycle_periods = cycle_times(2)-cycle_times(1);
     for i = 3:length(cycle_times)
         current_period = cycle_times(i) - cycle_times(i-1);
         cycle_periods = [cycle_periods current_period];
     end
% % 
% %     % Now we have cycle periods and since we start at 0 (or if we didn't start
% %     % at 0 we could offset by the amount of recording we do before zap. If we
% %     % record 2s before then 2+cycletime[0] = start time of second cycle
     offset = cycle_times(1);
     start_times = offset;
     for i = 1:length(cycle_periods)-1
         offset = offset + cycle_periods(i);
         start_times = [start_times offset];
     end
 
     frequencies = 1000./(cycle_periods);
%     % Now for each cycle period we can calculate the time from the start of the
%     % cycle to each successive spike that falls in that cycle duration
     spikephases = zeros(120, length(cycle_periods));
     spiketimes = zeros(120, length(cycle_periods));
     count = 0;
     flag = 0;
     cross_point = 0;
   
     for i = 1: length(start_times)-1
         start = find(t == start_times(i));
         finish = find(t == start_times(i+1));
         for j = start: finish-1
            % use these indicies to index vcut array, find the time of the spike
            % and normalize by the current cycle period i

            % only record the spike on the rising phase and ignore the falling
            % phase
            if(icut(j) > i_thr && flag == 1) 
                if (abs((t(j)-t(start)) - (cross_point-(t(j)-t(start)))) > 10)
                    count = count + 1; 
                    %store spike phase
                    spikephases(count, i) = (t(j)-t(start));
                    spiketimes(count, i) = (t(j)-t(start));
                end
                cross_point = t(j)-t(start);
                flag = 0;
            end
            if(icut(j) < i_thr && flag ==0)
                cross_point = t(j)-t(start);  
                flag = 1;
            end
         end
         spiketimes(:, i) = spikephases(:, i);
         spikephases(:, i) = spikephases(:, i)./cycle_periods(i);
         count = 0;
     end
% 
 
     avgisis = [];
     count = 1;
     num_isis = 0;
     isis = zeros(length(spiketimes(:, 1))-1, length(cycle_periods));
     % calculate ISIs for each cycle period
     for i = 1:length(cycle_periods)
         n = find(spiketimes(:, i) > 0)
         if(length(n) > 1)
             for j = 2:length(n)
                 isis(j-1, count) = spiketimes(j, i)-spiketimes(j-1, i);
             end
             % get average for all ISIs
             avgisis = [avgisis; mean(isis(1:n, count))];
             count = count + 1;
         else
            avgisis = [avgisis; NaN];
        end
     end
     
    spikephases(spikephases == 0) = NaN;
    spiketimes(spiketimes == 0) = NaN;
    
    avgif = 1./avgisis;
    figure(2)
    avgifsmooth = smooth(avgif,0.5,'loess');
    plot(frequencies, avgif, 'ob', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.99 .2 .23], 'MarkerSize', 5)
    hold on
    plot(frequencies, avgifsmooth, 'ob', 'MarkerEdgeColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.49 1 .63], 'MarkerSize', 5)
end

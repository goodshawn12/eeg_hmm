function [] = plotTimelock(axe, event, state_by_epoch, start_offset, Fs, win_len_sec, rt, rt_off, cmap, plot_thin)
    w = [1 1 1];
    k = [0 0 0];
    g = [0.8 0.8 0.8];
    thick = 2;
    thin = 0.01;

    if ismember(0, state_by_epoch)
        cmap = [1,1,1; cmap];
    end
    
    imagesc(state_by_epoch); hold on;
    if strcmp(event, '251')
        title = text('String','Time locked to lane-departure event');
        rt_x = rt - start_offset;
        rt_off_x = rt_off - start_offset;
        rt_linewidth = thick;
        rt_off_linewidth = thin;
        timelock_linecolor = k;
        rt_linecolor = w;
        rt_off_linecolor = g;
        if ~plot_thin
            plot(rt_x, 1:length(rt), 'linewidth', rt_linewidth, 'color', rt_linecolor);
        end
    elseif strcmp(event, '253')
        title = text('String','Time locked to car driver response start');
        rt_x = -start_offset - rt;
        rt_off_x = rt_off - start_offset;
        rt_linewidth = thick;
        rt_off_linewidth = thin;
        rt_linecolor = k;
        timelock_linecolor = w;
        rt_off_linecolor = g;
        if ~plot_thin
            plot(rt_x, 1:length(rt), 'linewidth', rt_linewidth, 'color', rt_linecolor);
        end
    elseif strcmp(event, '254')
        title = text('String','Time locked to driver response end');
        rt_x = -start_offset - rt_off;
        rt_off_x = -start_offset - rt_off + rt;
        rt_linewidth = thin;
        rt_off_linewidth = thin;
        rt_linecolor = k;
        rt_off_linecolor = w;
        timelock_linecolor = g;
    end
    xlabel = text('String','Time offset (sec)');
    ylabel = text('String', 'Epochs');
    
    x_ticks = 1:Fs:Fs*win_len_sec+1;
    x_tick_labels = cell(1, win_len_sec);
    x_tick_label_start = start_offset / Fs;    
    labels = x_tick_label_start:1:x_tick_label_start+win_len_sec;
    for i = 1:length(x_ticks)
        x_tick_labels{i} = num2str(labels(i));
    end

    set(axe,'Title', title);
    set(axe,'XLabel', xlabel);
    set(axe,'YLabel', ylabel);
    xticks(axe, x_ticks);
    xticklabels(axe, x_tick_labels);
    colormap(cmap);

    hold on,
    line([-start_offset, -start_offset], [1, length(rt)], 'linewidth', thick, 'color', timelock_linecolor);
    if plot_thin
        plot(rt_x, 1:length(rt), 'linewidth', rt_linewidth, 'color', rt_linecolor);
        plot(rt_off_x, 1:length(rt), 'linewidth', rt_off_linewidth, 'color', rt_off_linecolor);
    end
    
end
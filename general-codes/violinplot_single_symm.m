function h = violinplot_single_symm(data1, pos, width, color)
    if nargin < 2
        pos = 1;
        width = 0.3;
        color = [0, 0, 0];
    end

    % if size(data, 2) == 1
    %     data = repmat(data, 1, 2);
    % end

    if ~exist('data2')
        data2 = data1;
    end

    %% calculate kernel density estimate
    [density_1, value_1] = ksdensity(data1); %data(:, 1));
    % density_1 = density_1(value_1 >= min(data(:, 1)) & value_1 <= max(data(:, 1)));
    density_1 = density_1(value_1 >= min(data1) & value_1 <= max(data1));
    % value_1 = value_1(value_1 >= min(data(:, 1)) & value_1 <= max(data(:, 1)));
    value_1 = value_1(value_1 >= min(data1) & value_1 <= max(data1));
    % value_1(1) = min(data(:, 1));
    value_1(1) = min(data1);
    % value_1(end) = max(data(:, 1));
    value_1(end) = max(data1);

    % [density_2, value_2] = ksdensity(data(:, 2));
    [density_2, value_2] = ksdensity(data2);
    % density_2 = density_2(value_2 >= min(data(:, 2)) & value_2 <= max(data(:, 2)));
    density_2 = density_2(value_2 >= min(data2) & value_2 <= max(data2));
    % value_2 = value_2(value_2 >= min(data(:, 2)) & value_2 <= max(data(:, 2)));
    value_2 = value_2(value_2 >= min(data2) & value_2 <= max(data2));
    % value_2(1) = min(data(:, 2));
    value_2(1) = min(data2);
    % value_2(end) = max(data(:, 2));
    value_2(end) = max(data2);
    
    %% adjust width
    width = width/max([density_1, density_2]);
    
    % %% make violin plot
    % h = fill([pos+density_2*width pos-density_1(end:-1:1)*width], [value_2 value_1(end:-1:1)], color);
    
    % %% adjust violin properties
    % h.FaceAlpha = 0.25;
    % h.EdgeColor = [1, 1, 1];

    % ======================
    %% make violin plot
    % ======================
    % left
    h_l = fill([pos*ones(size(density_1)) pos-density_1(end:-1:1)*width], [value_1 value_1(end:-1:1)], color);
    h_l.FaceAlpha = 0.5;
    h_l.EdgeColor = [1, 1, 1];
    hold on
    h = 1;

    % right
    h_r = fill([pos+density_2*width pos*ones(size(value_2))], [value_2 value_2(end:-1:1)], color);
    h_r.FaceAlpha = 0.25;
    h_r.EdgeColor = [1, 1, 1];

    h = [h_l, h_r];

end
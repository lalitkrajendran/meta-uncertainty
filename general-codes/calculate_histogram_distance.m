function d = calculate_histogram_distance(y1, y2, bins, distance_method)

    % calculate histograms
    hist1 = histcounts(y1(:), bins, 'normalization', 'pdf');
    hist2 = histcounts(y2(:), bins, 'normalization', 'pdf');

    % calculate distance
    d = pdist2(hist1, hist2, str2func(distance_method));
end
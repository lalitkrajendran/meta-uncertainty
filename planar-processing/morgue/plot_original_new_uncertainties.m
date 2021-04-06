function plot_original_new_uncertainties(imx, imy, mcx, mcy, csx, csy, individual_method_array)
    figure
        
    subplot(1, 2, 1)
    plot(imx(1, :), imx(2, :), 'o')
    hold on
    plot(mcx(1, :), mcx(2, :), 'o')
    plot(csx(1, :), csx(2, :), 'o')
    plot([0 0.1], [0 0.1], 'color', [0 0 0 0.5])
    axis equal
    axis([0 0.1 0 0.1])
    box off
    xlabel('Original')
    ylabel('Sub')
    title('X')
    lgd = legend(individual_method_array, 'location', 'northoutside', 'orientation', 'horizontal');
    lgd.Position = [0.35 0.8 0.3732 0.0476];

    subplot(1, 2, 2)
    plot(imy(1, :), imy(2, :), 'o')
    hold on
    plot(mcy(1, :), mcy(2, :), 'o')
    plot(csy(1, :), csy(2, :), 'o')
    plot([0 0.1], [0 0.1], 'color', [0 0 0 0.5])
    axis equal
    axis([0 0.1 0 0.1])
    box off
    xlabel('Original')
    ylabel('Sub')
    title('Y')
    % legend(individual_method_array, 'location', 'northoutside', 'orientation', 'horizontal')
end
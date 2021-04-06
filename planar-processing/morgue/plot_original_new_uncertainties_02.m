function plot_original_new_uncertainties_02(imx, imy, mcx, mcy, csx, csy)
    figure
        
    subplot(1, 3, 1)
    plot([imx(1, :), imy(1, :)], [imx(2, :), imy(2, :)], 'o')
    hold on
    plot([0 0.2], [0 0.2], 'color', [0 0 0 0.5])
    axis equal
    axis([0 0.2 0 0.2])
    box off
    xlabel('Uncertainty, Original (pix.)')
    ylabel('Uncertainty, Sub (pix.)')
    title('IM')

    subplot(1, 3, 2)
    plot([mcx(1, :), mcy(1, :)], [mcx(2, :), mcy(2, :)], 'o')
    hold on
    plot([0 0.2], [0 0.2], 'color', [0 0 0 0.5])
    axis equal
    axis([0 0.2 0 0.2])
    box off
    xlabel('Uncertainty, Original (pix.)')
    ylabel('Uncertainty, Sub (pix.)')
    title('MC')

    subplot(1, 3, 3)
    plot([csx(1, :), csy(1, :)], [csx(2, :), csy(2, :)], 'o')
    hold on
    plot([0 0.2], [0 0.2], 'color', [0 0 0 0.5])
    axis equal
    axis([0 0.2 0 0.2])
    box off
    xlabel('Uncertainty, Original (pix.)')
    ylabel('Uncertainty, Sub (pix.)')
    title('CS')

    set(gcf, 'resize', 'off')
    set(gcf, 'position', [202   378   840   400])
    drawnow();

end
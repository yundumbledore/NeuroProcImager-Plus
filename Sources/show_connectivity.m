function show_connectivity(W_sequence)
    m1 = W_sequence(:,:,50); % show connectivity matrix when the subject was conscious
    m2 = W_sequence(:,:,800); % show connectivity matrix when the subject was unconscious

    minColorLimit = min([m1, m2],[],'all');
    maxColorLimit = max([m1, m2],[],'all');

    fig = figure('Position',[100 100 800 330]);
    subplot(1,2,1);
    imagesc(m1) 
    xlabel('Brain region index')
    ylabel('Brain region index')
    title('When the subject was conscious')

    subplot(1,2,2);
    imagesc(m2) 
    xlabel('Brain region index')
    ylabel('Brain region index')
    title('When the subject was unconscious')

    sgtitle('Inter-Regional Connectivity Estimates')

    h = axes(fig,'visible','off'); 
    c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]); % attach colorbar to h
    c.Title.String = 'Strength'; % add colorbar title
    clim(h,[minColorLimit,maxColorLimit]);             % set colorbar limits
end
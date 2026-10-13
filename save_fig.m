function save_fig(fig, name, cfg)
% SAVE_FIG  Export a figure as a 300-dpi PNG into cfg.figDir (README 1.2).
%   save_fig(gcf, 'fig_pca_trajectory', cfg)
    exportgraphics(fig, fullfile(cfg.figDir, [name '.png']), 'Resolution', 300);
end

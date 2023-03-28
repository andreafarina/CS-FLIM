function [] = FormatFig(hfig,fname)
picturewidth = 5.7*0.8; % 4.2 set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',8*0.8) % 6 adjust fontsize to your document
set(gca, 'FontName', 'Arial')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
end
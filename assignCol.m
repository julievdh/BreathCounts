% assignCol: Assign colours to species

c = colormap(viridis);
cind = round(linspace(1,length(viridis),length(unique(spp(~isnan(spp)))))); % choose incides for colors based on how many spp
call = c(cind,:); % colours per species entry

% sort colours by weight of spp automatically
[~,ind] = sort(sp_wt);
%cind = cind(ind);

for fl = 1:length(files)
    if ~isnan(files(fl).spp) == 1
        if     strfind(files(fl).tag,'tt') == 1 % tursiops
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'hp') == 1 % phocoena
            files(fl).col = c(cind(ind(files(fl).spp)),:);
        end
        if     strfind(files(fl).tag,'pw') == 1 % short finned pilot whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'gm') == 1 % long finned pilot whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'sw') == 1 % sperm whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'zc') == 1 % cuvier's beaked whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'md') == 1 % baird's beaked whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'ha') == 1 % northern bottlenose whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'dl') == 1 % beluga
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'bb') == 1 % antarctic minke whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'mn') == 1 % humpback whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'eg') == 1 % right whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'bp') == 1 % fin whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'bw') == 1 % blue whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'gg') == 1 % risso's
            files(fl).col = c(cind(files(fl).spp),:);
        end
        if     strfind(files(fl).tag,'er') == 1 % gray whale
            files(fl).col = c(cind(files(fl).spp),:);
        end
    else files(fl).col = [0 0 0];
    end
end


%% check colours and weights
figure(99), clf, hold on
for i = 1:length(sp_ind)
    hold on
    bar(i,sp_wt(i),'Facecolor',files(sp_ind(i)).col)
    text(i-0.2,30,files(sp_ind(i)).tag(1:2),'color','w')
end
set(gca,'yscale','log','xtick',[1:16])

hold off

for i = 1:length(sp_ind)
    bar(phyorder(i),sp_wt(i),'Facecolor',files(sp_ind(i)).col)
    hold on
    text(phyorder(i)-0.2,30,files(sp_ind(i)).tag(1:2),'color','w')
    text(phyorder(i)-0.2,100,num2str(phyorder(i)),'color','k')
end
set(gca,'yscale','log','xtick',[1:16])

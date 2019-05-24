% assignCol: Assign colours to species

c = colormap(viridis);
cind = round(linspace(1,length(viridis),length(unique([files(:).spp])))); % choose incides for colors based on how many spp
call = c(cind,:); % colours per species entry

for fl = 1:length(files)
    if     strfind(files(fl).tag,'tt') == 1 % tursiops
        files(fl).col = c(cind(files(fl).spp),:);
    end
    if     strfind(files(fl).tag,'hp') == 1 % phocoena
        files(fl).col = c(cind(files(fl).spp),:);
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
    if     strfind(files(fl).tag,'ba') == 1 % minke whale
        files(fl).col = c(cind(files(fl).spp),:);
    end
    if     strfind(files(fl).tag,'mn') == 1 % humpback whale
        files(fl).col = c(cind(files(fl).spp),:);
    end
    %if     strfind(files(fl).tag,'eg') == 1 % right whale
    %    files(fl).col = c(cind(files(fl).spp),:);
    %end
    if     strfind(files(fl).tag,'bp') == 1 % fin whale
        files(fl).col = c(cind(files(fl).spp),:);
    end
    if     strfind(files(fl).tag,'bw') == 1 % blue whale
        files(fl).col = c(cind(files(fl).spp),:);
    end
    
end
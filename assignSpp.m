% assignspp - assign species code

for fl = 1:length(files)
    if strfind(files(fl).tag,'tt') == 1 % tursiops
        files(fl).spp = 2;
    end
    if     strfind(files(fl).tag,'hp') == 1 % phocoena
        files(fl).spp = 1;
    end
    if     strfind(files(fl).tag,'pw') == 1 % short finned pilot whale
        files(fl).spp = 5;
    end
    if     strfind(files(fl).tag,'gm') == 1 % long finned pilot whale
        files(fl).spp = 7;
    end
    if     strfind(files(fl).tag,'sw') == 1 % sperm whale
        files(fl).spp = 10;
    end
    if     strfind(files(fl).tag,'zc') == 1 % cuvier's beaked whale
        files(fl).spp = 8;
    end
    if     strfind(files(fl).tag,'md') == 1 % baird's beaked whale
        files(fl).spp = 6;
    end
    if     strfind(files(fl).tag,'ha') == 1 % Northern bottlenose whale
        files(fl).spp = 9;
    end
    if     strfind(files(fl).tag,'dl') == 1 % beluga whale
        files(fl).spp = 4;
    end
    if     strfind(files(fl).tag,'ba') == 1 % minke whale
        files(fl).spp = 11;
    end
    if     strfind(files(fl).tag,'gg') == 1 % risso's dolphin
        files(fl).spp = 3;
    end
    if     strfind(files(fl).tag,'mn') == 1 % humpback whale
        files(fl).spp = 14;
    end
    if     strfind(files(fl).tag,'er') == 1 % gray whale
        files(fl).spp = 13;
    end
    if     strfind(files(fl).tag,'eg') == 1 % right whale
        files(fl).spp = 12;
    end
    if     strfind(files(fl).tag,'bp') == 1 % fin whale
        files(fl).spp = 15;
    end
    if     strfind(files(fl).tag,'bw') == 1 % blue whale
        files(fl).spp = 16;
    end
    if     strfind(files(fl).tag,'bb') == 1 % antarctic minke
        files(fl).spp = 11;
    end
end

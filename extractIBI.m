% extract info from file structure 

for i = 1:length(files)
    if isempty(files(i).tdiff) == 0
        % calculate duration of resp audit in hours
        files(i).dur = (files(i).resp(end) - files(i).resp(1))/3600;
        if files(i).dur >= 6
            mnIBI(i) = mean(files(i).tdiff); % mean IBI
            stdIBI(i) = std(files(i).tdiff); % STD
            mdIBI(i) = median(files(i).tdiff); % median IBI
            modeIBI(i) = mode(files(i).tdiff); % mode
            
            % TLC(i) = 0.135*(files(i).wt).^0.92; % estimate from Kooyman
            % if isfield(files(i),'swim_ct') == 1
            mnf(i) = 60/mnIBI(i); % mean frequency
            mdf(i) = 60/mdIBI(i); % median frequency
            % end
        else
            mnIBI(i) = NaN; % mean IBI
            stdIBI(i) = NaN; % STD
            mdIBI(i) = NaN; % median IBI
            modeIBI(i) = NaN; % mode
            files(i).dur = NaN; % nan these out so they're not used in the calcs below
            files(i).wt = NaN; 
            files(i).spp = NaN; 
        end
        else files(i).dur = NaN;
            files(i).spp = NaN; 
    end
end

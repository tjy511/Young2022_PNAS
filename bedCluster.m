function bedDepths = bedCluster(bedsday, clustNum)

% Find total ranges of possible bed depth clusters
try
    % Restrict number of clusters to 3 most probablee number 
    [A,~] = kmeans(bedsday,clustNum-1);
    [B,~] = kmeans(bedsday,clustNum);
    [C,~] = kmeans(bedsday,clustNum+1);
    
    % Obtain means and variances for each cluster
    for ii = 1:clustNum-1
        Amean(ii) = mean(bedsday(A==ii));
        Avar(ii) = var(bedsday(A==ii));
    end
    
    for ii = 1:clustNum
        Bmean(ii) = mean(bedsday(B==ii));
        Bvar(ii) = var(bedsday(B==ii));
    end
    
    for ii = 1:clustNum+1
        Cmean(ii) = mean(bedsday(C==ii));
        Cvar(ii) = var(bedsday(C==ii));
    end
    
    % Identify most likely number of clusters through observing variances.
    % Returns the mean of the identified cluster. 
    thresh = 0.001;
    if all(Avar<thresh) == 1
        bedDepths = [Amean NaN NaN];
    elseif all(Bvar<thresh) == 1
        bedDepths = [Bmean NaN];
    else
        bedDepths = Cmean;
    end
catch
    % Returns NaN if number of clusters were incorrectly identified. 
    bedDepths = nan(1,clustNum);
end
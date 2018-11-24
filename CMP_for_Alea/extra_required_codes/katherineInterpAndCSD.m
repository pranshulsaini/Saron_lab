% This script is for interpolation of channels in Katherine's data.
function katherineInterpAndCSD()
    % data
    dataset = '/Volumes/EEG STORAGE/EEG_Analyses/CPT_R2/pre/eeg_analysis/prestim_alpha/id38_epoch_alpha_target3.set';
    newsfp = '/Volumes/EEG STORAGE/EEG_Analyses/CPT_R2/pre/eeg_analysis/prestim_alpha/pre_id38_sfp.sfp';
    EEG = pop_loadset('filename', dataset);
    EEG.chanlocs = readlocs(newsfp);
    badCh = {'Z1'};
    EEGbad = [];    
	oldsfp = '/Volumes/EEG STORAGE/EEG_Analyses/CPT_R2/pre/preSPid38_pre2_day1_09_08_07_dig.sfp';
    EEGbad.chanlocs = readlocs(oldsfp);
    EEGbad.chanlocs = EEGbad.chanlocs(match_str({EEGbad.chanlocs.labels},badCh));
    for i = 1:1:length(EEGbad.chanlocs)
       if EEGbad.chanlocs(i).labels(1) == 'X'
           EEGbad.chanlocs(i).labels(1) = 'A';
       elseif EEGbad.chanlocs(i).labels(1) == 'Y'
           EEGbad.chanlocs(i).labels(1) = 'B';
       elseif EEGbad.chanlocs(i).labels(1) == 'Z'
           EEGbad.chanlocs(i).labels(1) = 'C';
       end
    end
    interpEEG = [];
    data = [];
    for epoch = 1:1:size(EEG.data,3)
        EEGtemp = EEG;
        EEGtemp.data = EEG.data(:,:,epoch);
        newEEG = spline_interpolation(EEGtemp,EEGbad);
        if epoch == 1
            interpEEG = newEEG;
        end
        data(:,:,epoch) = newEEG.data;
    end
    interpEEG.data = data;
    
    interpAndCsdEEG = ApplyCSDTransform(interpEEG); 
    
end
function csdEEG = ApplyCSDTransform(interpEEG)
    csdEEG = interpEEG;
    el = interpEEG.chanlocs;

    phi = [el.sph_phi]';
    theta = [el.theta]';
    
    phiT = 90 - phi;                    % calculate phi from top of sphere
    theta2 = (2 * pi * theta) / 360;    % convert degrees to radians
    phi2 = (2 * pi * phiT) / 360;
    [x,y] = pol2cart(theta2,phi2);      % get plane coordinates
    xy = [x y];
    xy = xy/max(max(xy));               % set maximum to unit length
    xy = xy/2 + 0.5;                    % adjust to range 0-1

    mel.lab = {el.labels}';
    
    mel.phi = phi;
    mel.theta = theta;
    mel.xy = xy;

    [G,H] = GetGH(mel); 
    lambda = 1e-6;
    
    for tr = 1:1:size(csdEEG.data,3)
        csdEEG.trial{tr} = data2.trial{tr}-mean(data2.trial{tr},2)*ones(1,size(data2.trial{tr},2));
        data2.trial{tr} = CSD(data2.trial{tr},G,H,lambda);    
    end
    
end
function newEEG = spline_interpolation(EEG, EEGbad)
    % get theta, rad of electrodes
    % ----------------------------
    xelec = [ EEG.chanlocs.X ];
    yelec = [ EEG.chanlocs.Y ];
    zelec = [ EEG.chanlocs.Z ];
    rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
    xelec = xelec./rad;
    yelec = yelec./rad;
    zelec = zelec./rad;
    xbad = [ EEGbad.chanlocs.X ];
    ybad = [ EEGbad.chanlocs.Y ];
    zbad = [ EEGbad.chanlocs.Z ];
    rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
    xbad = xbad./rad;
    ybad = ybad./rad;
    zbad = zbad./rad;
    [tmp1 tmp2 tmp3 badchansdata] = spheric_spline_int( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data);
    EEGbad.data = badchansdata;
    
    newEEG=EEG;
    newEEG.chanlocs(end+1:end+size(badchansdata,1)) = EEGbad.chanlocs;
    newEEG.data(end+1:end+size(badchansdata,1),:) = EEGbad.data;
    newEEG.nbchan = newEEG.nbchan + length(EEGbad.chanlocs);
end
function [xbad, ybad, zbad, allres] = spheric_spline_int( xelec, yelec, zelec, xbad, ybad, zbad, values)

    newchans = length(xbad);
    numpoints = size(values,2);

    Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
    Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);

    % compute solution for parameters C
    % ---------------------------------
    meanvalues = mean(values); 
    values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero

    values = [values;zeros(1,numpoints)];
    C = pinv([Gelec;ones(1,length(Gelec))]) * values;
    clear values;
    allres = zeros(newchans, numpoints);

    % apply results
    % -------------
    for j = 1:size(Gsph,1)
        allres(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));        
    end
    allres = allres + repmat(meanvalues, [size(allres,1) 1]);
end
% compute G function
function g = computeg(x,y,z,xelec,yelec,zelec)

    unitmat = ones(length(x(:)),length(xelec));
    EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                    (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                    (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

    g = zeros(length(x(:)),length(xelec));
    %dsafds
    m = 4; % 3 is linear, 4 is best according to Perrin's curve
    for n = 1:7
        L = legendre(n,EI);
        g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
    end
    g = g/(4*pi);    
end
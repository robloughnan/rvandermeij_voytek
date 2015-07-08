function lnspectra = rmr_findlinespectra(dat,fsample,searchrange,param)

% This function iteratively finds line spectra in electrophysiological
% recordings, with the purpose of filtering these line spectra afterwards
% in a different pipeline. Function outputs are peak frequencies of the line
% spectra combined with their filter half-bandwidths that suppress the line
% spectra successfully. Specification of the filters used for this purpose 
% are given in the output as well (and can given as input), as well as
% the appropriate length of the edges that should be cut off after
% filtering (in the case of IIR, an estimate based on 98% of the abs of the
% approximated impulse response). 
% 
% Filters are meant to applied to the entire recording regardless of on which 
% channel line spectra occur. Because of this, peaks and bandwidths are 
% determined using all channels. This is done by averaging their processed
% PSDs.
% 
% Spectral analysis performed to detect peaks is done using a Hanning 
% tapered Welch window. The resulting PSD (log10) is than filtered itself to
% remove any slow trends. The resulting processed PSD is then Z-valued
% (mean/std calculed for all channels simultaneously) and averaged over channels.
% The result processed PSD is then used to find peaks.
% 
% The peak finding algorithm consists of two iterative loops, having the following logic:
% 1) get PSD in normal/log space
% 2) filter PSD to flatten it
% 3) z-value and average over channels
% 4) threshold and determine peaks as center of contiguous frequency bins above threshold
% 5) if peaks exist:
%   a) filter data
%   b) get PSD in normal/log space
%   c) filter PSD to flatten it
%   d) z-value and average over channels (using z-value from 3)
%   e) if peaks are not fully suppressed, increase bandwidth and return to a
% 6) return to 1 for next pass using filtered data from last step 
% 
% 
% 
% 
%                Input:
%                         dat = Nchan X Nsample spatio-temporal matrix 
%                     fsample = scalar, sampling rate in Hz
%                 searchrange = 1x2 vector, frequency range in Hz to search in [freqmin freqmax]
%                       param = structure, containing additional inputs
%               param.zthresh = scalar, Z-value to use as threshold for catching peaks (default = 2)
%              param.welchwin = scalar, length in seconds of Welch window to use for PSD (default = 10)
%             param.bandwstep = scalar, half bandwidth in Hz, stepsize for determining filter bandwidth (default = 0.25)
%                 param.taper = string, taper to use for spectral estimation (default = 'hanning')
%              param.filttype = string, type of bandstop filter  'but' (default),  'fir', etc (FieldTrip-style name)
%               param.filtord = scalar, order of bandstop filter (default = 2)
%               param.filtdir = string, direction of bandstop filter 'twopass' (default), 'onepass', etc (FieldTrip-style name)
%
%
% 
%                Output:
%              lnspectra.peak = 1xN vector of peak frequencies in Hz
%         lnspectra.halfbandw = 1xN vector of half-bandwidths used for filtering (peak +/- halfbandw)
%          lnspectra.filttype = string, filter type of bandstop filter used
%           lnspectra.filtord = scalar, order of  bandstop filter used
%           lnspectra.filtdir = string, filter direction of bandstop filter used
%        lnspectra.edgeartlen = scalar, length in seconds of the edge artifact to be removed
%           lnspectra.origpsd = 1xN vector, PSD before filtering, for inspection purposes
%           lnspectra.filtpsd = 1xN vector, PSD after filtering, for inspection purposes
% 
%
%  (this function depends on FieldTrip, https://github.com/fieldtrip)
%  
%

%
% Copyright (C) 2015, Roemer van der Meij, roemervandermeij AT gmail DOT com
%

% set defaults 
if ~isfield(param,   'zthresh'),     param.zthresh     = 2;            end
if ~isfield(param,   'welchwin'),    param.welchwin    = 10;           end
if ~isfield(param,   'bandwstep'),   param.bandwstep   = 0.25;         end
if ~isfield(param,   'taper'),       param.taper       = 'hanning';    end
if ~isfield(param,   'filttype'),    param.filttype    = 'but';        end
if ~isfield(param,   'filtord'),     param.filtord     = 2;            end
if ~isfield(param,   'filtdir'),     param.filtdir     = 'twopass';    end

% sanity checks
if size(dat,2)<size(dat,1)
  error('dat has more channels than samples')
end


%%%%%%%%%%%%%%%%%%%%%
% iteratively find peaks, and iteratively adjust bandwidth of the peaks found
peaksremaining = true;
peaks     = [];
bandwidth = [];
pspecpeaks = [];
pspecbandw = [];
pspecprpow = [];
count = 0;
while peaksremaining
  count = count + 1;
    
  
  % get pow
  [pow, freq] = getpow(dat,fsample,searchrange,param.welchwin,param.taper);
  % process pow
  [procpow,zparam] = processpow(pow,freq,[],peaks,bandwidth);
  
  % save origpow for output if first pass and save round specific processed pow
  if count == 1 
    origpow = pow;
  end
  pspecprpow{count} = procpow;
  
  % find  peaks
  [peaks, bandwidth] = findpeaks(procpow,freq,param.zthresh,param.bandwstep);
  if ~isempty(peaks)
    disp(['found ' num2str(numel(peaks)) ' line spectra in pass ' num2str(count)])
    for ipeak = 1:numel(peaks)
      disp(['line spectra at ' num2str(peaks(ipeak)) 'Hz +/- ' num2str(bandwidth(ipeak)) 'Hz'])
    end
  else
    % no more peaks found, quit
    peaksremaining = false;
    disp(['no more peaks found in pass ' num2str(count)])
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%
  % iteratively increase bandwidth till all peaks are gone
  peakgone = false(size(peaks));
  while ~all(peakgone)
    
    % filter data
    filtdat = dat;
    for ipeak = 1:numel(peaks)
      % apply a bandstop filter
      disp(['applying filter for peak at ' num2str(peaks(ipeak)) 'Hz +/- ' num2str(bandwidth(ipeak)) 'Hz'])
      filtdat = ft_preproc_bandstopfilter(filtdat, fsample, [peaks(ipeak)-bandwidth(ipeak) peaks(ipeak)+bandwidth(ipeak)], param.filtord, param.filttype, param.filtdir);
    end
    
    % get pow and process it, using same zval-ling as used initially
    [pow, freq] = getpow(filtdat,fsample,searchrange,param.welchwin,param.taper);
    procpow = processpow(pow,freq,zparam);
    
    % assess residual peak presence per peak, report in peaksgone, and increase bandwidth if necessary
    for ipeak = 1:numel(peaks)
      currpeak  = peaks(ipeak);
      currbandw = bandwidth(ipeak);
      currind   = find(freq>=(currpeak-currbandw) & freq<=(currpeak+currbandw));
      if any(procpow(currind)>param.zthresh)
        changed = true;
        bandwidth(ipeak) = bandwidth(ipeak) + param.bandwstep;
        disp(['increased bandwidth for line spectra at ' num2str(peaks(ipeak)) 'Hz to +/-' num2str(bandwidth(ipeak)) 'Hz'])
      else
        peakgone(ipeak) = true;
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%
  
  % save pass-specific peaks and bandwidth
  pspecpeaks{count} = peaks;
  pspecbandw{count} = bandwidth;
  % set filtered dat to dat for next round
  dat = filtdat;

end
%%%%%%%%%%%%%%%%%%%%%
% gather peaks and bandwidths
peaks     = cat(2,pspecpeaks{:});
bandwidth = cat(2,pspecbandw{:});


% determine approximate length of impulse response at 98% of impulse response (determined by sum of abs)
edgeartlen = NaN(1,numel(peaks));
for ipeak = 1:numel(peaks)
  fakedat = zeros(1,round(10*fsample));
  fakedat(round(5*fsample)) = 1;
  impresp  = ft_preproc_bandstopfilter(fakedat, fsample, [peaks(ipeak)-bandwidth(ipeak) peaks(ipeak)+bandwidth(ipeak)], param.filtord, param.filttype, param.filtdir);
  impresp  = impresp ./ sum(abs(impresp));
  edgeartlen(ipeak) = ((find(cumsum(abs(impresp))>0.98,1)-round(5*fsample))*2)./fsample;
end
edgeartlen = max(edgeartlen);


% create output structure
lnspectra = [];
lnspectra.peak       = peaks;
lnspectra.halfbandw  = bandwidth;
lnspectra.filttype   = param.filttype;
lnspectra.filtord    = param.filtord;
lnspectra.filtdir    = param.filtdir;
lnspectra.edgeartlen = edgeartlen;
lnspectra.origpsd    = origpow;
lnspectra.filtpsd    = pow; % latest pow is maximally filtered pow




%%%%%%%%%%%%%%%%%
% show some plots of results and process
% 1) original PSD in log/norm space with found peaks highlighted and filtered spectrum on top of it
figure('numbertitle','off','name','original and filtered mean PSD in log/normal space')
hold on
l1 = plot(freq,mean(log10(lnspectra.origpsd)),'color',rgb('green'));
for ipeak = 1:numel(peaks)
  begfreq = peaks(ipeak)-bandwidth(ipeak);
  endfreq = peaks(ipeak)+bandwidth(ipeak);
  ind = find(freq>=begfreq & freq<=endfreq);
  l2 = plot(freq(ind),mean(log10(lnspectra.origpsd(:,ind))),'color',rgb('red'));
end
l3 = plot(freq,mean(log10(lnspectra.filtpsd)),'color',rgb('blue'));
xlabel('frequency (Hz)')
ylabel('log mean (over channels) power')
legend([l1 l2 l3],'orignal PSD','identified as line spectrum','filtered PSD');

% 2) z-valued the original PSD in log/norm space with found peaks highlighted and filtered spectrum on top of it
figure('numbertitle','off','name','processed PSD of each pass in log/normal space')
npass = numel(pspecpeaks);
for ipass = 1:npass
  subplot(1,npass,ipass)
  hold on
  currpeaks = pspecpeaks{ipass};
  currbandw = pspecbandw{ipass};
  currprpow = pspecprpow{ipass};
  l1 = plot(freq,currprpow,'color',rgb('blue'));
  for ipeak = 1:numel(currpeaks)
    begfreq = currpeaks(ipeak)-currbandw(ipeak);
    endfreq = currpeaks(ipeak)+currbandw(ipeak);
    ind = find(freq>=begfreq & freq<=endfreq);
    l2 = plot(freq(ind),currprpow(ind),'color',rgb('red'));
  end
  l3 = line([freq(1) freq(end)],[param.zthresh param.zthresh],'linestyle','--','color',rgb('dark green'));
  xlabel('frequency (Hz)')
  ylabel('log mean (over channels) power')
  legend([l1 l2 l3],'processed PSD','identified as line spectrum','threshold');
  title(['PSD and identified peaks of pass ' num2str(ipass)])
  ylim = get(gca,'ylim');
  set(gca,'ylim',[-3 max([ylim(2) 5])])
  set(gca,'xlim',[freq(1) freq(end)])
end
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% SUB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate power
function [pow, freq] = getpow(dat,fsample,searchrange,welchwin,taper)

% set basics
nsample = size(dat,2);
nchan   = size(dat,1);

%%%% get PSD using a welch approach WILL BE SUBFUNCTION LATER
% determine the welch windows to loop over, 75% overlap
nsampwelch   = round(welchwin .* fsample); % samples in a welch window
welchwinstep = round(nsampwelch ./ 4); % steps in between windows
welchindbeg  = 1:welchwinstep:nsample;
welchindend  = nsampwelch:welchwinstep:nsample;
welchindbeg  = welchindbeg(1:numel(welchindend));
welchind  = [welchindbeg; welchindend].';
nwelchwin = size(welchind,1);
% determine freq axis and freq bins to use for searching
freqboilim = round([searchrange(1) searchrange(2)] ./ (fsample ./ nsampwelch)) + 1;
freqboi    = freqboilim(1):1:freqboilim(2);
freq       = (freqboi-1) ./ (welchwin);
% get PSD
pow = zeros(nchan,numel(freqboi));
for iwelch = 1:nwelchwin
  
  % cut out and taper
  tap = feval(taper,nsampwelch);
  tap = tap ./ sqrt(sum(tap.^2));
  currind = welchind(iwelch,:);
  currdat = dat(:,currind(1):currind(2));
  currdat = bsxfun(@times,currdat,tap.');
  
  % get power
  fftdat = fft(currdat,[], 2);
  currpow = abs(fftdat(:,freqboi)).^2; % select bins
  
  % calc running mean
  pow = pow + (currpow ./ nwelchwin);
end
% the above is equivalent to pwelch(dat.',nsampwelch,welchwinstep*3,freq,fsample).', but pwelch used in this way is horribly slow
% the former also allows for more control (tapering)
%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Log power, fit slope, remove it from non log space, z-value them all together, and mean
function [procpow, zparam] = processpow(pow,freq,varargin)

% get variable input arg
if nargin == 2
  zparam = [];
elseif nargin == 3
    zparam = varargin{1};
elseif nargin == 5
  zparam    = varargin{1};
  peaks     = varargin{2};
  bandwidth = varargin{3};
else
  error('internal error: nargin should 2, 3, or 5')
end

% log pow and freq axis
logpow  = log10(pow);
logfreq = log10(freq);

% filter it with a very low highpass. The annoying trick I applied here, the variant of mirror padding, is because I couldn't 
% get a decent filter response at low (approximate) order without inducing edge artifacts
% create a fake time using the frequency resolution
fakefsample = 1./(mean(diff(freq))); % to make things easy, consider 1Hz as 1s
faketimelen = (freq(end)-freq(1));
% make filter which killes everything whose frequeny doens't fits less than X times in the whole fake time signal
% there's an cumbersome trade-off here: the smaller X, the more wiggles stay in there, and fitting capability is reduced.
% but the higher X, the more the 'slower' edges of broad line spectra get pushed down, which also affects fitting capability
hpfreq  = 6/faketimelen;
% filter settings are hardcoded and not meant to be changed
filttype = 'but';
filtord  = 2;
filtdir  = 'twopass';
% apply filtering with variant of 'mirror' padding
padlen  = 100*fakefsample; % 100Hz 
prepad  = logpow(:,1:padlen)-repmat(logpow(:,padlen)-logpow(:,1),[1 padlen]); % prepad is start of powspctrm, copied, and realigned to the start of the spectrum. This is to maintain the slope in the start of the spectrum, which we want this pad out smoothly to aid filtering
postpad = logpow(:,end:-1:end-padlen); % postpad is a mirror of the edge
filtpow = [prepad logpow postpad];
filtpow = ft_preproc_highpassfilter(filtpow, fakefsample, hpfreq, filtord, filttype, filtdir);
filtpow = filtpow(:,(padlen+1):end-(padlen+1)); % remove padding

% keeping here for playing around, not helping currently
% % fit using rmr_robustfit
% for ichan = 1:nchan
%   offschi = rmr_robustfit(logfreq,filtpow(ichan,:));
%   filtpow(ichan,:) = filtpow(ichan,:) - (offschi(1) + offschi(2).*logfreq);
% end

% zval it, and output mean/std as zparam
if isempty(zparam)
  % calculate zparams on parts of the spectrum that are not peaks (which shoot up/down before/after bandstop filtering)
  tmpfiltpow = filtpow;
  if ~isempty(peaks)
    for ipeak = 1:numel(peaks)
      currpeak  = peaks(ipeak);
      currbandw = bandwidth(ipeak);
      currind   = find(freq>=(currpeak-currbandw) & freq<=(currpeak+currbandw)); % robust to boundaries
      tmpfiltpow(:,currind) = NaN;
    end
  end
  meanlpow = nanmean(filtpow(:));
  stdlpow  = nanstd(filtpow(:));
  zparam   = [meanlpow stdlpow];
end
procpow  = (filtpow - zparam(1)) ./ zparam(2);
procpow  = mean(procpow,1);
%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find peaks and determine bandwidth
function [peaks, bandwidth] = findpeaks(procpow,freq,zthresh,startbandw)

% find threshold crossings
passthresh = procpow>zthresh;
% find contigious points and select center of these
if any(passthresh)
  passind = [];
  passind(:,1) = find(diff([0,passthresh,0]) == 1);
  passind(:,2) = find(diff([0,passthresh,0]) == -1) - 1;
  peaks = freq(round(mean(passind,2)));
  % get bandwidth start values
  bandwidth = freq(passind(:,2))-freq(passind(:,1));
  bandwidth(bandwidth<startbandw) = startbandw;
  % ensure maximal accuracy of .1Hz
  peaks     = round(peaks*100)/100;
  bandwidth = round(bandwidth*100)/100;
else
  peaks     = [];
  bandwidth = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%









function playground



[poworg, freqorg] = getpow(orgdat,fsample,searchrange,10,param.taper);
[poworg075, freq075] = getpow(orgdat,fsample,searchrange,0.075,param.taper);
[filtpow, freqorg] = getpow(dat,fsample,searchrange,10,param.taper);
[filtpow075, freq075] = getpow(dat,fsample,searchrange,.075,param.taper);

figure
hold on
plot(freqorg,mean(log10(poworg)),'color',rgb('light blue'))
plot(freq075,mean(log10(poworg075)),'color',rgb('blue'))
plot(freqorg,mean(log10(filtpow)),'color',rgb('light green'))
plot(freq075,mean(log10(filtpow075)),'color',rgb('green'))
legend('original PSD 10s Welch','original PSD 0.075s Welch','filtered PSD 10s Welch','filtered PSD 0.075s Welch')
































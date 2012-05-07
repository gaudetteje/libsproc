function fftgroup(varargin)
% FFTGROUP  Calculates the frequency spectrum for a large group of time
% series

% Author:   Jason Gaudette
% Company:  Naval Undersea Warfare Center (Newport, RI)
% Phone:    401.832.6601
% Email:    gaudetteje@npt.nuwc.navy.mil
% Date:     20060928
%

Navg = 4;   % number of pulses to average together
winname = 'blackman7';
winstart = 0.001;   % ignore first 1 ms
winlength = 0.006;  % use a 6 ms window
tmin = 1;   % minimum number of seconds in time series to look at

tic;    % start timer

% get directory to search
if (nargin > 0)
    wdir = char(varargin(1));
else
    wdir = pwd;
end
if (wdir(end) ~= filesep)
    wdir = [wdir filesep];
end

% log output to file
diary([wdir sprintf('fftgroup_%s.log', datestr(now,'yyyymmdd'))]);
disp(datestr(now))

% compile search results
flist = findfiles(wdir, '\.dlm$');
nfiles = length(flist);
if ~nfiles
    error('No files found in directory "%s"', wdir);
    diary 'off';
end

disp(sprintf('Found %d files...', nfiles));

% initialize figure window for plotting
%close all hidden;
%fh = figure('Visible', 'off');% , 'DockControls', 'off');

for fnum=1:nfiles
    % clear variables for next iteration
    clear ts fd
    
    % parse string for file and path names
    fname = char(flist(fnum));
    ind = max(strfind(fname,filesep));
    fdir = fname(1:ind);
    fname_dlm = fname(ind+1:end);
    fname_mat = [fname_dlm(1:end-4) '.mat'];
    fname_fft = [fname_dlm(1:end-4) ' FFT ' winname '.mat'];
    fname_fig = [fname_fft(1:end-4) '.fig'];
    
    % check for previously calculated FFT
    if exist([fdir fname_fft], 'file')
        disp(sprintf('[%d] Frequency spectrum already calculated!  Bypassing "%s"', fnum, [fdir fname_mat]));
        continue;
    end
    
    % load DLM file
    if ~exist([fdir fname_dlm], 'file')
        disp(sprintf('[%d] No delimited pulse location file found!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    try
        %pulses = load([fdir fname_dlm], '-ascii');
        pulses = load([fdir fname_dlm]);
    catch
	pulses = [];
    end

    if ~length(pulses)
        disp(sprintf('[%d] No pulses specified in DLM file!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    
    % load MAT file
    load([fdir fname_mat]);
    if exist('fs', 'var')
        disp(sprintf('[%d] Frequency spectrum already calculated!  Bypassing "%s"', fnum, [fdir fname_mat]));
        continue;
    end

    if ~exist('ts', 'var')
        disp(sprintf('[%d] No time series data struct found!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    
    % store default settings for later reference
    fd.fname = ts.fname;
    fd.fs = ts.fs;
    fd.navg = Navg;
    fd.winname = winname;
    fd.winstart = winstart;
    fd.winlength = winlength;
    fd.tmin = tmin;
    fd.points = round(fd.winlength*fd.fs);
    
    pmin = find(ts.time(pulses) > tmin, 1)-1;    % find the first pulse to start with
    
    % verify we have enough pulses to work with
    if ~length(pmin) || (length(pulses) < (pmin + fd.navg))
        disp(sprintf('[%d] Not enough pulses available in time series!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    fprintf('[%d] Processing "%s" using %s window\n', fnum, fname_dlm, winname);
    
    % take Navg distinct pulses
    for N = 1:Navg
        fd.ind(N,:) = [1:fd.points] + round(fd.winstart*ts.fs) + pulses(pmin+N);
        fd.data(N,:) = fftpulse(ts.data(fd.ind(N,:)), 'linear', fd.winname)';
    end
    
    % calculate frequency sequence (in units of kHz)
    fd.freq = [-fd.points/2 : fd.points/2-1] * fd.fs / (fd.points * 1000);
    %fd.freq = linspace(0,1,fd.points/2+1) * (fd.fs/2000);

    
    % calculate average over Navg FFTs
    fd.avg = mean(fd.data, 1);
    fd.avgDB = 20*log10(fd.avg);
    fd.units = 'dBV';
    
    % plot data to figure
    plot(fd.freq, fd.avgDB);
    title([sprintf('%d point %s FFT: ', fd.points, fd.winname) fd.fname(1:end-4)]);
    xlabel('Frequency (kHz)');
    ylabel(fd.units);
%     ymin = floor(1.05 * min(fd.avgDB)/20)*20;
%     ymax = ceil(1.05 * max(fd.avgDB)/20)*20;
    axis([0 fd.freq(end) -140 0]); grid on;
    
%%% WIN32
    % save plot of averaged spectrum
    try
        set(fh, 'Visible', 'on');
        saveas(fh, [fdir fname_fig], 'fig');
        set(fh, 'Visible', 'off');
    catch
        warning('Could not save figure to file: "%s"\n', fname_fig);
        continue
    end
    
    % save spectrum data to mat file
    try
        save([fdir fname_fft], 'fd', '-MAT');
    catch
        warning('Could not save data to file: "%s"\n', fname_mat);
        continue
    end
%%% LINUX
    % save plot of averaged spectrum
%    try
%        print([fdir fname_eps], '-dpsc2');
%	if ~exist([fdir fname_eps], 'file')
%		pause(0.2)
%	end
%	system(sprintf('convert %s -rotate 90 %s', [fdir fname_eps], [fdir fname_fig]));
%    catch
%        warning('Could not save figure to file: "%s"\n', fname_fig);
%    end
    
    % append spectrum data to mat file
%    try
%        save('-mat', [fdir fname_mat], 'ts', 'fd');
%    catch
%        warning('Could not save data to file: "%s"\n', fname_mat);
%        continue
%    end

end

fprintf('Finished processing %d files in %.0f seconds (%.2f minutes).\n\n', nfiles, toc, toc/60);
diary off;

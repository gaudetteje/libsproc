function locgroup(varargin)
% LOCGROUP  Locates CW pulses for a large group of time series data

% Author:   Jason Gaudette
% Company:  Naval Undersea Warfare Center (Newport, RI)
% Phone:    401.832.6601
% Email:    gaudetteje@npt.nuwc.navy.mil
% Date:     20060928
%

gamma = 200;

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
diary([wdir sprintf('locgroup_%s.log', datestr(now,'yyyymmdd'))]);
disp(datestr(now))

% compile search results
flist = findfiles(wdir, '\.mat$');
nfiles = length(flist);
if ~nfiles
    error('No files found in directory "%s"', wdir);
    diary 'off';
end

disp(sprintf('Found %d files...', nfiles));

for fnum=1:nfiles
    clear ts pulses
    
    % parse string for file and path names
    fname_mat = char(flist(fnum));
    ind = max(strfind(fname_mat,filesep));
    fdir = fname_mat(1:ind);
    fname_mat = fname_mat(ind+1:end);
    fname_dlm = [fname_mat(1:end-4) '.dlm'];
    
    % verify file does not already exist
    if exist(fname_dlm, 'file')
        disp(sprintf('[%d] File already processed!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    
    % load MAT file into the workspace
    load([fdir fname_mat]);
    if ~exist('ts', 'var')
        disp(sprintf('[%d] No time series data struct found!  Bypassing "%s"', fnum, fname_mat));
        continue;
    end
    
    % get frequency of TX pulse for matched filter in locpulse
%     if exists(ts.fc, 'var')
%         fc = ts.fc;
%     else
        ind = regexpi(ts.fname, '_');
        ts.fc = str2num(ts.fname(ind(1)+1:ind(2)-2)) * 1000;
        if ~length(ts.fc)
            % estimate fc by maximum FFT level rounded to nearest 5kHz in range
            %ts.fc = 0;
            continue;
        end
%     end
    
    % locate CW pulses and store in struct
    pulses = locpulse(ts.data, ts.fs, ts.fc, 0.001, 0.3);
    disp(sprintf('[%d] Located %d pulses in file "%s"', fnum, length(pulses), ts.fname));
    
    % save pulse locations in DLM file
    try
        %save([fdir fname_dlm], 'pulses', '-ASCII', '-DOUBLE');
        save('-ascii', [fdir fname_dlm], 'pulses');
    catch
        warning('Could not save data to file: "%s"\n', fname_dlm);
    end
end

disp(sprintf('Finished processing %d files in %.0f seconds (%.2f minutes).\n\n', nfiles, toc, toc/60));
diary off;


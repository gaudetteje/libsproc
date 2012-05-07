function Z = ambiguity_function(fs,fc,signal_1,signal_2)

mode = 2; % mode = 1 makes (rows,cols) = (time,Doppler), 2 transposes this to (Doppler,time)
plot_results = 1; % 1 to plot image of ambiguity diagram, 0 to not
scale_results = 0; %1 to normalize maximum correlation to 0 dB FOR PLOTTING ONLY
%
% If only one signal is supplied, use autocorrelation
%
try
  signal_2;
catch
  signal_2 = signal_1;
end
signal_1 = signal_1(:);
signal_2 = signal_2(:);
%
% If the user gave us a single value for fc, compute the Doppler frequency
% shifts from this value.
%
% If the user gave us a (n x 2) or (2 x n) matrix the assume it holds
% Doppler (knots) in one row or column and frequuency (Hz) in the other row
% or column
%
if max(size(fc)) == 1
  %
  % Compute Doppler frequency shifts based on fc
  %
  v_increment = 0.2; % knots
  max_v = 40; % knots
  c = 1500; % meters/sec
  c_knots = c * 3600/1852; % convert c from meters/sec to knots (nautical miles per hour)

  v = [-fliplr(v_increment:v_increment:max_v) 0:v_increment:max_v];
  frequencies = ((1 + v/c_knots) ./ (1 - v/c_knots))*fc - fc;
else
  %
  % Use user-supplied Doppler frequency shifts
  %
  if min(size(fc)) ~= 2
    error('AMBIGUITY_FUNCTION: fc must either be a scalar or, at most, an (n x 2) or (2 x n) matrix');
  else
    %
    % Ensure one column is Doppler and the other frequency
    %
    [rows,cols] = size(fc);
    if rows == cols
      error('AMBIGUITY_FUNCTION: Must have at least 3 Doppler frequencies.')
    elseif cols > rows
      fc = fc';
    end
    v = fc(:,1);
    frequencies = fc(:,2);
  end
end

if length(signal_1) > length(signal_2)
  t = (0:length(signal_1)-1)/fs;
  time_delay = (-length(signal_1)+1:1:length(signal_1)-1)/fs;
else
  t = (0:length(signal_2)-1)/fs;
  time_delay = (-length(signal_2)+1:1:length(signal_2)-1)/fs;
end
t = t(:);
Z = zeros(2*max(length(signal_1),length(signal_2))-1,length(frequencies));

if length(signal_1) > length(signal_2)
  for freq_count = 1:length(frequencies)
    Z(:,freq_count) = xcorr(signal_2,(signal_1) .* exp(j*2*pi*frequencies(freq_count)*t));
  end
else
  for freq_count = 1:length(frequencies)
    Z(:,freq_count) = xcorr(signal_1,(signal_2) .* exp(j*2*pi*frequencies(freq_count)*t));
  end
end
%
% Format the output to handle the users desired (time,Doppler) orientation
%
if mode == 1
  % do nothing, just here for argumant checking purposes
elseif mode == 2
  Z = Z.';
else
  error('AMBIGUITY_FUNCTION: MODE must be 1 or 2!')
end


if plot_results == 1
  if mode == 1
    if scale_results == 1
      imagesc(v,time_delay*1000,20*log10(abs(Z./(max(max(Z))))));
    else
      imagesc(v,time_delay*1000,20*log10(abs(Z)));
    end
    caxis([-60 0])
    xlabel('Velocity (knots)')
    ylabel('Time Delay (ms)')
  else
    if scale_results == 1
      imagesc(time_delay*1000,v,20*log10(abs(Z./(max(max(Z))))));
    else
      imagesc(time_delay*1000,v,20*log10(abs(Z)));
    end
    caxis([-60 0])
    xlabel('Time Delay (ms)')
    ylabel('Velocity (knots)')
  end
end

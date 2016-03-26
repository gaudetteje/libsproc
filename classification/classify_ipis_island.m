function varargout = classify_ipis_island(x)
% CLASSIFY_IPIS_ISLAND  use island and stability criteria to label each
% pulse in a sonar sound group
%
% [y1, y2, y3, y4] = classify_ipis_island(IPI)  takes an array of 
% inter-pulse intervals in IPI and returns the index to the first pulse in
% each singlet, doublet, etc. according to the island and stability
% criteria of Kothari, et al.
%
% classify_ipis_island(IPI)  plots the results when no output parameters
% are specified.
%
%
% To Do:
% - Adaptively extend analyis beyond 4 consecutive pulses (should identify
% landing/terminal buzz pulse trains as one event).  Currently these are
% picked out as singlets.
% - Use input parameters to change stability threshold (i.e., alpha, gamma)
%
%
% Reference:
%     Kothari, N. B., Wohlgemuth, M. J., Hulgard, K., Surlykke, A., & Moss, 
%     C. F. (2014). Timing matters: sonar call groups facilitate target
%     localization in bats. Frontiers in Physiology, 5, 168?13.
%     http://doi.org/10.3389/fphys.2014.00168

% Author:   Jason Gaudette
% Company:  Naval Undersea Warfare Center (Newport, RI)
% Phone:    401.832.6601
% Email:    jason.e.gaudette@navy.mil
% Date:     20160326

% pick stability threshold for strobe groups > 2
alpha = 1.2;
gamma = 0.05;

% data check
x = x(:)';   % force into row vector
assert(all(x >= 0), 'Event times must be positive!');

% find and remove last pulse events (ID by NaN)
idxLast = isnan(x);
x(idxLast) = 1e3;       % force post-IPI to larger than 100 ms

% define event times from IPIs
t = [0 cumsum(x)] .* 1e-3;

% define pre-post sequences
x_n1 = [Inf x(1:end-1)];
x_p1 = [x(2:end) inf];
x_p2 = [x(3:end) inf(1,2)];
x_p3 = [x(4:end) inf(1,3)];

mu1 = (x + x_p1)./2;
mu2 = (x + x_p1 + x_p2)./3;

% find doublets
idx2 = find((x*alpha <= x_n1) & (x*alpha <= x_p1));

% find triplets
idx3a = find((x*alpha <= x_n1) & (x*alpha <= x_p2));
idx3b = find((abs(mu1 - x)./mu1 < gamma) & (abs(mu1 - x_p1)./mu1 < gamma));
idx3 = intersect(idx3a, idx3b);

% find quadruplets
idx4a = find((x*alpha <= x_n1) & (x*alpha <= x_p3));
idx4b = find((abs(mu2 - x)./mu2 < gamma) & (abs(mu2 - x_p2)./mu2 < gamma));
idx4 = intersect(idx4a, idx4b);

% infer singlets from remaining pulses
idx1 = unique([idx2 idx2+1 idx3 idx3+1 idx3+2 idx4 idx4+1 idx4+2 idx4+3]);
idx1 = setxor(idx1, (1:numel(x)));

%% plot events on timeline
if nargout == 0
    figure
    plot(t,zeros(size(t)),'.k')
    hold on
    plot(t(idx1),zeros(size(idx1)),'ob')
    plot(t(idx2),zeros(size(idx2)),'og')
    plot(t(idx3),zeros(size(idx3)),'or')
    plot(t(idx4),zeros(size(idx4)),'om')

    % plot IPIs over events
    plot(t(1:end-1),x,'*k')

    xlabel('Time (s)')
    ylabel('post-IPI (ms)')
    title('Classification of Pulse Types by IPI')

    classtypes = {'post-IPI','singlet','doublet','triplet','quadruplet'};
    if isempty(idx4), classtypes = classtypes(1:end-1); end
    if isempty(idx3), classtypes = classtypes(1:end-1); end
    if isempty(idx2), classtypes = classtypes(1:end-1); end
    legend(classtypes,'IPI event values')
    
    % ensure no duplicate labels
    assert(all(size(unique([idx2 idx3 idx4])) == size([idx2 idx3 idx4])), 'duplicate sound group labels found');

    % ensure every data point is labeled
    assert(numel(x) == numel(idx1) + numel(idx2)*2 + numel(idx3)*3 + numel(idx4)*4, 'missing labels found')
end

% otherwise return indices
for n=1:nargout
    eval(sprintf('varargout{%d} = idx%d;',n,n));
end

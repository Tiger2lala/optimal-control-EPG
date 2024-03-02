function b1map = b1convert(Adj, refV, ROI)
% convert B1 map from Adj file to epg code format
% mz407

B1abs = reshape(Adj.S, Adj.image_m, Adj.image_n, Adj.coils, Adj.slices);
B1abs = permute(B1abs, [3 1 2 4]);

if nargin < 3
    B1list(:,:) = B1abs(:, Adj.W);
else
    B1list(:,:) = B1abs(:, ROI);
end

% So far it's in uT/V
% Need to turn into % of prescribed FA
% FA(deg) * 500Hz * 1ms / 180deg / gamma / (0.5 * Vref * 1ms) = B1 (T/V)
gamma = 42.577e6;
B1fadeg = B1list * 0.5 * refV * gamma * 180 / 500 /1e6;

b1perc = B1fadeg / 90;

% pre divide by sqrt(8); add CP phase
b1map = b1perc / sqrt(8) .* exp(-[0:7]/8*2*pi*1j).';
b1map = b1map.';




clear all;
close all;



%dat = importdata('data/caporaso_m3_raw.dat'); %microbiome (M3)
dat = importdata('data/plankton_bacteria_raw.dat'); %plankton (bacteria)

abn = dat.data; %table
sps = dat.textdata(2:end); %species

time = abn(1, :); %time
abn = abn(2:end, :); %abundances

%normalization
[sp tm] = size(abn);
abn = abn ./ repmat(sum(abn), sp, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delete single peak time series
j = [];
for i = 1:length(abn(:, 1))
	clear L num g
	[L num] = bwlabel(abn(i, :));
	for ii = 1:num
		g(ii) = length(find(L == ii));
	end

	if length(find(g > 1)) > 0
		j = [j i];
	end
end

abn = abn(j, :);
[sp tm] = size(abn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I(t)
It = zeros(sp, tm - 1);
for i = 1:sp
	aux = [];
	for t = 2:tm
		aux = [aux (abn(i, t) - abn(i, t - 1))/std(abn(i, :))];
	end
	It(i, :) = aux;
end

[sp tm] = size(It);
rIt = reshape(It, sp*tm, 1);
h = histogram(rIt, 100);
vals = h.Values;
vals = vals/max(vals);
medios = [];
for i = 2:length(h.BinEdges)
	medios = [medios mean([h.BinEdges(i-1), h.BinEdges(i)])];
end

figure(1)
plot(medios, log10(vals), 'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ccdf
[sp tm] = size(abn);
rabn = reshape(abn, sp*tm, 1);
rabng0 = rabn(find(rabn > 0));

[xbin ybin] = ccdf(rabng0, 2500);

figure(2)
plot(log10(xbin), log10(ybin), '-xk')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% taylor
mx = [];
sx = [];
for i = 1:sp
	mx = [mx mean(abn(i, :))];
	sx = [sx std(abn(i, :))];
end

figure(3)
loglog(mx, sx, 'x')




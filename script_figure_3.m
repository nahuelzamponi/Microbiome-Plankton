


close all
clear all



dat = importdata('data/caporaso_m3_raw.dat'); %microbiome (M3)
%dat = importdata('data/plankton_bacteria_raw.dat'); %plankton (bacteria)

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

% avalanches/events
[sp, t] = size(abn);
aval = zeros(sp, t);
for i = 1:sp
	ab = abn(i, :);
	[L, num] = bwlabel(ab);

	clear s ini
	for ii = 1:num
		x = find(L == ii);
		s(ii) = trapz(time(x), ab(x), 2);
		ini(ii) = min(x);
	end

	ze = zeros(1, t);
	ze(ini) = s;
	aval(i, :) = ze;
end

auxi = [];
xaux = [];
kk = 1;
for i = 1:length(aval(:, 1))
	clear ava

	if sum(aval(i, :)) > 0
		ava = aval(i, :);
		auxi = [auxi ava];
		xaux = [xaux time];
	end
end

auxi = auxi(find(auxi > 0));
xaux = xaux(find(auxi > 0));

figure(1)
plot(xaux, log10(auxi), '*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ccdf
bin = 0.0001;
au = bin;
while au <= 0.05
	bin = [bin bin(length(bin))*2];
	au = bin(length(bin));
end

figure(2)
ava = reshape(aval, sp*t, 1);
pos = find(ava > 0);
av = ava(pos);
[bins H] = ccdf(av, 5000000, 1);
loglog(bins, H, '--*')
hold on
for i = 1:length(bin)
	xline(bin(i))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stationary
dt = 10; %wdw size
fr = {};
for k = 1:length(bin)
	ini = 1;
	fin = dt;

	center = [];
	aux = [];
	cont = 1;
	while fin < t
		center = [center (time(ini) + time(fin))/2];
		LL = reshape(aval(:, ini:fin), sp*dt, 1);
		aux(cont) = length(find(LL >= bin(k)));
		cont = cont + 1;
	
		ini = fin + 1;
		fin = fin + dt;
	end

	fr{k} = aux;
end

figure(3)
for i = 1:length(fr)
	plot(center, log10(fr{i}), '-x')
	hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% histograms

figure(4)
b = {};
h = {};
for i = 1:length(fr)
	hi = histfit(fr{i}, 10, 'kernel')
	b{i} = hi(2).XData;
	h{i} = hi(2).YData;
end

for i = 1:length(b)
	clear bb bbb hh hhh
	bb = b{i}(find(b{i} > 0));
	hh = h{i}(find(b{i} > 0));

	bbb = bb(find(hh >= 1));
	hhh = hh(find(hh >= 1));

	plot(log10(bbb), log10(hhh), '-')
	hold on
end
%xlim([-1 2]);
%ylim([0 0.8]);

figure(5)
normo = [];
b = {};
h = {};
for i = 1:length(fr)
	hi = histfit(fr{i}, 10, 'kernel');
	b{i} = hi(2).XData;
	h{i} = hi(2).YData;
end

for i = 1:length(b)
	clear bb bbb hh hhh posi
	bb = b{i}(find(b{i} > 0));
	hh = h{i}(find(b{i} > 0));

	bbb = bb(find(hh >= 1));
	hhh = hh(find(hh >= 1));
	
	[val posi] = max(hhh);
	
	normo(i) = bbb(posi);

	bbb = bbb/bbb(posi);
	hhh = hhh/max(hhh);

	plot(log10(bbb), log10(hhh), '-')
	hold on
end
%xlim([-1.5 0.5]);
%ylim([-0.7 0.1]);




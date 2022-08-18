


close all
clear all



dat = importdata('data/caporaso_m3_genera.dat'); %microbiome (M3)

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

N = [];
max_eigen = [];

sorteo = {};
rank = {};

rnk = 5:-1:1; % number of taxa included (phylum, class, order, family, genus)
for ii = 1:length(rnk)
	taxa = [];
	for i = 1:length(sps)
		clear aux
		aux = strsplit(sps{i}, '_');

		auxi = [];
		for kk = 1:rnk(ii)
			auxi = [auxi string(aux{kk})];
		end
		taxa = [taxa join(auxi, '_')];
	end

	utax = unique(taxa);
	umat = zeros(length(utax), length(time));
	for i = 1:length(utax)
		j = find(taxa == utax(i));
		umat(i, :) = sum(abn(j, :), 1);
	end

	[sp tm] = size(umat);
	N(ii) = sp;

	%eigen
	r = corr(umat');
	ee = real(eig(r));
	ee = ee(find(ee > 10^-7)); % discard extremely low values
	max_eigen(ii) = max(ee);
	sorteo{ii} = sort(ee, 'descend');
	rank{ii} = [1:length(sorteo{ii})]';
end

figure(1)
for i = 1:length(sorteo)
	loglog(rank{i}, sorteo{i}, '--o')
	hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collapse
tx1 = [];
for i = 1:length(rank)
	tx1 = [tx1 max(rank{i})];
end
scaling_exp = [0.66274 0.63623 0.58119 0.69259 0.99963]; 

figure(2)
for i = 1:length(rank)
	loglog(rank{i}./tx1(i)^0.7074, sorteo{i}.*(rank{i}./tx1(i)).^scaling_exp(i), '-o')
	hold on;

	nrank{i} = rank{i}./tx1(i)^0.7074;
	nsorteo{i} = sorteo{i}.*(rank{i}./tx1(i)).^scaling_exp(i);
end

grid on
axis([0.01 25 0.000001 50]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inset
ejex = log10(N);
ejey = log10(max_eigen);
p1 = polyfit(ejex, ejey, 1);
f = polyval(p1, ejex);
p1

figure(3)
plot(log10(N), log10(max_eigen), 'o', ejex, f, '--r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Renormalization group (meshulam et al)
N = [];
max_eigen = [];

sorteo = {};
rank = {};
nrank = {};

dats = abn;
[sp t] = size(dats);
Xk{1} = dats;
r = corr(dats');
ee = real(eig(r));
ee = ee(find(ee > 10^-7));
max_eigen(1) = max(ee);
sorteo{1} = sort(ee, 'descend');
rank{1} = [1:length(sorteo{1})]';
nrank{1} = rank{1}/max(rank{1});
N(1) = sp;

for i = 2:6
	xk = zeros(floor(sp/2), t);
	r = corr(dats');
	rr = r;
	rr(:, :) = 1;
	rr = tril(rr)*-1;
	rr = rr + triu(r, 1);

	ki = 0;
	while ki < floor(sp/2)
		maximum = max(max(rr));
		[ii, jj] = find(rr == maximum);

		xij = dats(ii(1), :) + dats(jj(1), :);
		nxij = xij/mean(xij(find(xij > 0)));

		ki = ki + 1;
		xk(ki, :) = nxij;

		rr(ii(1), :) = -1;
		rr(jj(1), :) = -1;
		rr(:, ii(1)) = -1;
		rr(:, jj(1)) = -1;
	end

	r = corr(xk');
	ee = real(eig(r));
	ee = ee(find(ee > 10^-7));
	max_eigen(i) = max(ee);
	sorteo{i} = sort(ee, 'descend');
	rank{i} = [1:length(sorteo{i})]';
	nrank{i} = rank{i}/max(rank{i});

	Xk{i} = xk;
	dats = xk;
	[sp t] = size(dats);
	N(i) = sp;
end

figure(4)
for i = 1:length(sorteo)
	loglog(rank{i}, sorteo{i}, '--o')
	hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%collapse
tx2 = [];
for i = 1:length(rank)
	tx2 = [tx2 max(rank{i})];
end
scaling_exp = [0.65131 0.66606 0.60807 0.68498 0.91041 0.66917];

figure(5)
for i = 1:length(rank)
	loglog(rank{i}.*tx2(i).^-0.7439, sorteo{i}.*(rank{i}./tx2(i)).^scaling_exp(i), '-o') %use the same mu for everything
	hold on;

	nrank{i} = rank{i}.*tx2(i)^-0.7439;
	nsorteo{i} = sorteo{i}.*(rank{i}./tx2(i)).^scaling_exp(i);
end

grid on
axis([0.001 50 0.000001 5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inset
ejex = log10(N);
ejey = log10(max_eigen);
p1 = polyfit(ejex, ejey, 1);
f = polyval(p1, ejex);
p1

figure(6)
plot(log10(N), log10(max_eigen), 'o', ejex, f, '--r')














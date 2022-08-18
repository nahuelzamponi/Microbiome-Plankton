


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

th1 = 0;
th2 = 0;

T = [];
S = [];
aval = {};
t_ava = {};
cont = 1;
for i = 1:sp
	ab = abn(i, :);
	abb = ab;
	abb(find(abb > th1)) = 1; %binarisation
	[L num] = bwlabel(abb);

	clear t s
	for ii = 1:num
		x = find(L == ii);
		t(ii) = length(time(min(x)):time(max(x)));
		s(ii) = trapz(time(x), ab(x), 2);

		if length(x) >= th2
			aval{cont} = ab(x);
			t_ava{cont} = time(x);
			cont = cont + 1;
		end
	end

	T = [T t];
	S = [S s];
end

%remove s = 0 from S and T
Tp = T(find(S > 0));
Sp = S(find(S > 0));



xn = [];
xnp = [];
for i = 1:length(aval)
	xn = [xn aval{i}(1, 1:(length(aval{i}) - 1))];
	xnp = [xnp aval{i}(1, 2:length(aval{i}))];
end

xx = 0.00001:0.001:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th = 0.0001;

%activity
S = [];
for i = 1:sp
	ab = abn(i, :);
	abb = ab;
	abb(find(abb >= th)) = 1; %binarisation
	[L num] = bwlabel(abb);

	clear s
	for ii = 1:num
		clear x
		x = find(L == ii);
		s(ii) = trapz(time(x), ab(x), 2);
	end

	S = [S s];
end

%remove s = 0 from S and T
St = S(find(S > 0));

figure(2)
[bins H] = ccdf(St, 0.00001, 1);
loglog(bins, H, '--*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%rest
R = [];
for i = 1:sp
	ab = mat(i, :);
	abb = ab;
	abb(find(abb >= th)) = 2; %binarisation
	abb(find(abb < th)) = 1; %binarisation
	abb(find(abb == 2)) = 0; %binarisation
	[L num] = bwlabel(abb);

	if sum(abb) > 0
		clear s
		for ii = 1:num
			clear x
			x = time(find(L == ii));
			s(ii) = length(min(x):max(x));
		end

		R = [R s];
	end
end

%remove s = 0 from S and T
Rth = R(find(R > 0));

figure(3)
[bins H] = ccdf(Rth, 0.00001, 1);
loglog(bins, H, '--*')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0 = mean(xnp(find(xn == min(xn))));

base = 2;

maxbin = round(logb(max(xn), base));
un = sort(unique(xnp));
minbin = round(logb(un(2), base)); %un(1) is used as v0
minbin = minbin - 2;

k = sort(xn);

i = minbin;
bin = [];
ii = [];
ff = [];
while i <= maxbin
	ii = base^i;
	ff = base^(i + 1);
	bin = [bin ff];
	i = i + 1;
end

m = [xn(:) xnp(:)];

[~, idx] = sort(m(:, 1));
sort_m = m(idx, :); 

clear mxnp sxnp
mxnp = mean( sort_m( find(sort_m(:, 1) <= bin(1)), 2 ) );
sxnp = std( sort_m( find(sort_m(:, 1) <= bin(1)), 2 ) );
for k = 2:length(bin)
	mxnp = [mxnp mean(sort_m(find((sort_m(:, 1) <= bin(k)) & (sort_m(:, 1) > bin(k - 1))), 2))];
	sxnp = [sxnp std(sort_m(find((sort_m(:, 1) <= bin(k)) & (sort_m(:, 1) > bin(k - 1))), 2))];
end

figure(1)
loglog(xn, xnp, '.')
hold on
loglog(xx, xx, '--r')
loglog(bin, mxnp, '-og')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

av_max = length(Sp);
delta = 0.0005;
special_delta = 0.001;

aux = v0;

l = 0;
tt = 2;
while l < av_max
	clear vtp pos

	if aux(tt - 1) == v0
		vtp = sort_m(find((sort_m(:, 1) <= (aux(tt - 1) + special_delta)) & (sort_m(:, 1) > (aux(tt - 1) - special_delta))), 2);
		pos = randi(length(vtp));
		aux = [aux vtp(pos)];
	else
		vtp = sort_m(find((sort_m(:, 1) <= (aux(tt - 1) + delta)) & (sort_m(:, 1) > (aux(tt - 1) - delta))), 2);
		pos = randi(length(vtp));
		aux = [aux vtp(pos)];
	end

	if vtp(pos) < th
		pausa = Rth(randi(length(Rth)));
		aux = [aux zeros(1, pausa)];
		l = l + 1;
	end

	tt = tt + 1;
end
syn_abs = aux;



time = 0:(length(syn_abs) - 1);
S = [];

ab = syn_abs;
abb = ab;
abb(find(abb > 0)) = 1; %binarisation
[L num] = bwlabel(abb);

clear s
for ii = 1:num
	x = find(L == ii);
	s(ii) = trapz(time(x), ab(x), 2);
end

S = [S s];

%remove s = 0 from S and T
Sth = S(find(S > 0));

figure(4)
clear bins H
[bins H] = ccdf(Sth, 0.00001, 1);
loglog(bins, H, '--*')



















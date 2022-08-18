


clear all;
close all;



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

% avalanche duration and size
[sp tm] = size(abn);
th1 = 0; %offset
th2 = 0; %av duration

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



subplot(2, 2, 1);
[t_bin t_cnt] = ccdf(Tp, 1000000, 0);
%[t_bin t_cnt] = ccdf(Tp, 500, 1);

ejex = log10(t_bin(2:25));
ejey = log10(t_cnt(2:25));
p1 = polyfit(ejex, ejey, 1);
p1
f = polyval(p1, ejex);
plot(log10(t_bin), log10(t_cnt), '-o', ejex, f, '-r')
xlabel('T')
ylabel('P(T)')



subplot(2, 2, 2);
[s_bin s_cnt] = ccdf(Sp, 1000000, 0);
%[s_bin s_cnt] = ccdf(Sp, 1000, 1);

ejex = log10(s_bin(10:100));
ejey = log10(s_cnt(10:100));
%ejex = log10(s_bin(2:100));
%ejey = log10(s_cnt(2:100));
p1 = polyfit(ejex, ejey, 1);
p1
f = polyval(p1, ejex);
plot(log10(s_bin), log10(s_cnt), '-o', ejex, f, '-r')
%xlim([0 4])
xlabel('S')
ylabel('P(S)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scaling T v <S>
TS = zeros(length(Tp(1, :)), 2);
TS(:, 1) = Tp;
TS(:, 2) = Sp;
TS = sortrows(TS, 1);

base = 1.5;
minbin = fix(logb(min(Tp), base));
maxbin = round(logb(max(Tp), base));

i = minbin;
delta = 0.5;

bin = minbin;
while i <= maxbin
	ii = base^i;
	ff = base^(i + delta);
	bin = [bin round(ff, 2)];
	i = i + delta;
end

mS = mean(TS(find(TS(:, 1) >= minbin & TS(:, 1) < bin(1)), 2));
sdS = std(TS(find(TS(:, 1) >= minbin & TS(:, 1) < bin(1)), 2));
for i = 2:length(bin)
	mS = [mS mean(TS(find(TS(:, 1) >= bin(i - 1) & TS(:, 1) < bin(i)), 2))];
	sdS = [sdS std(TS(find(TS(:, 1) >= bin(i - 1) & TS(:, 1) < bin(i)), 2))];
end

subplot(2, 2, 3);
bin = bin(find(mS > 0));
mS = mS(find(mS > 0));

ejex = log10(bin);
ejey = log10(mS);

p1 = polyfit(ejex, ejey, 1);
p1
f = polyval(p1, ejex);
plot(log10(bin), log10(mS), 'o', ejex, f, '--r')
xlabel('T')
ylabel('<S>(T)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shape collapse
paso = 0.01;
bins = 0:paso:1;
mnX = [];
sdnX = [];

th_min = 5;
th_max = 101;
aval2 = {};
t_ava2 = {};
%tt = [];
k = 1;
for i = 1:length(aval)
	if (length(aval{i}) >= th_min) & (length(aval{i}) < th_max)
		aval2{k} = [0 aval{i}/mean(aval{i}) 0];
		t_ava2{k} = 0:1/(length(aval{i}) + 1):1;
		%tt(k) = max(t_ava{i}) - min(t_ava{i});
		k = k + 1;
	end
end

nX = [];
tT = [];
for i = 1:length(aval2)
	%nX = [nX aval2{i}]; %(tt(i)^0.85)];
	nX = [nX aval2{i}];
	tT = [tT t_ava2{i}];
end

mnX(1) = mean(nX(find(tT == 0)));
sdnX(1) = std(nX(find(tT == 0)));
for i = 2:length(bins)
	mnX(i) = mean(nX(find(tT > bins(i - 1) & tT <= bins(i))));
	sdnX(i) = std(nX(find(tT > bins(i - 1) & tT <= bins(i))));
end



subplot(2, 2, 4);
%plot(tT, nX, 'LineWidth', 0.1, 'Color',[0.75 0.75 0.75])
%hold on
plot(bins, mnX, 'o')%, 'LineWidth', 1.5)
ylim([0.5 1.5])
xlabel('t/T')
%ylabel('<S>/T^m')
ylabel('x(t/T)/<x(t/T)>')




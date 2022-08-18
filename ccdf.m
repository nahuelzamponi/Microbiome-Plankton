function [bns frqs] = ccdf(dat, bins, cum);
%input cum = 1 if you want your instagram cumulative

if ~exist('cum', 'var')
	% third parameter does not exist, so default it to something
	cum = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute and normalise histogram
[vals edgs] = histcounts(dat, bins);

fr = vals;
bi = edgs;
bb = [];
kk = 1;
for i = 2:length(bi)
	bb(kk) = (bi(i) + bi(i - 1))/2;
	kk = kk + 1;
end

if cum == 0
	clear bns freq
	bns = bb(find(fr > 0));
	freq = fr(find(fr > 0));
	frqs = freq/max(freq);
else %cum = 1
	bns = bb;
	freq1 = fr/sum(fr);
	%flip the array
	freq2 = fliplr(freq1);
	%make the cumulative sum starting from the largest value
	freq3 = cumsum(freq2);
	%flip back the array to plot starting from zero
	frqs = fliplr(freq3);
end

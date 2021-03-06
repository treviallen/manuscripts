% NBCC2015_RTGM.m

% Parse files

locs02 = 'NBCC2015Loc_mean_hazcurves_02.csv';
locs10 = 'NBCC2015Loc_mean_hazcurves_10.csv';

locsfile = locs02

probs = [0.02	0.01375	0.01	0.00445	0.0021	0.001	0.0005	0.000404	0.0002	0.0001];
[lon lat par s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 ref loc prov] = ...
     textread(locsfile,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s','delimiter',',','headerlines',4);
     
% concat SAs

SAs = [s1 s2 s3 s4 s5 s6 s7 s8 s9 s10];

% get RTGMs

RTGMs = [];
RiskCoeffs = [];
Iplot = 0;
HazardCurve.AFEs = probs';
for i=1:length(loc)
	if ref(i) == 40 % for testing
		HazardCurve.SAs = SAs(i,:)'	
		disp(['Calculating RTGM for ',loc(i)]);
		[ RTGM, RiskCoefficient ] = RTGM_Calculator_Ver131017( HazardCurve, Iplot );
		RTGMs = [RTGMs RTGM];
		RiskCoeffs = [RiskCoeffs RiskCoefficient];	
	end
end

% export RTGM values
header = 'REF,LON,LAT,RTGM,RISKCOEFF'

data = [ref lon lat RTGMs' RiskCoeffs']

if locsfile == locs02
	outfile = 'NBCC2015Loc_RTGM_02.csv';
else
	outfile = 'NBCC2015Loc_RTGM_10.csv';
end
dlmwrite(outfile, header, 'delimiter','');
dlmwrite(outfile, data, 'delimiter',',', '-append');


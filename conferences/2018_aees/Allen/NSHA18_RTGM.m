% NSHA18_RTGM.m

% Parse files

locs02 = 'interp_haz_curves_SA(0.2).csv';
locs10 = 'interp_haz_curves_SA(1.0).csv';
locsPGA = 'interp_haz_curves_PGA.csv';

locsfile = locs10;

probs = [0.02	0.01375	0.01	0.00445	0.0021	0.001	0.0005	0.000404	0.0002	0.0001];
[lon lat s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 places] = ...
     textread(locsfile,'%f%f%f%f%f%f%f%f%f%f%f%f%s','delimiter',',','headerlines',1);
     
% concat SAs

SAs = [s1 s2 s3 s4 s5 s6 s7 s8 s9 s10];

% get RTGMs

RTGMs = [];
RiskCoeffs = [];
Iplot = 0;
HazardCurve.AFEs = probs';
for i=1:length(places)
	%if strcmp(places{i},'Adelaide') > 0 % for testing
		HazardCurve.SAs = SAs(i,:)'	;
		disp(['Calculating RTGM for ',places{i}]);
		[ RTGM, RiskCoefficient FragilityCurves ] = RTGM_Calculator_Ver131017( HazardCurve, Iplot );
		RTGMs = [RTGMs RTGM];
		RiskCoeffs = [RiskCoeffs RiskCoefficient ];	
	%end
end

% export RTGM values
header = ['PLACE,LON,LAT,RTGM,RISKCOEFF' char(10)];

%places = char(places);
%data = [places lon lat RTGMs' RiskCoeffs'];

if strcmp(locsfile, locs02) > 0
	outfile = 'NSHA18_RTGM_SA(0.2).csv';
elseif strcmp(locsfile, locs10) > 0
	outfile = 'NSHA18_RTGM_SA(1.0).csv';
else
	outfile = 'NSHA18_RTGM_PGA.csv';
end

% make out txt
outtxt = header;
for i=1:length(places)
    line = [places{i},',',num2str(lon(i)),',',num2str(lat(i)),',', ...
            num2str(RTGMs(i)),',',num2str(RiskCoeffs(i)),char(10)];
    outtxt = [outtxt line];
end
    
dlmwrite(outfile, outtxt, 'delimiter','');
%dlmwrite(outfile, data, 'delimiter',',', '-append');








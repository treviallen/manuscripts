function [] = Plot_RTGM_Calculation_Iterations( HazardCurve, AFE4UHGM, UHGM, ...
                  RTGMi, FragilityCurves, FRAGILITY_AT_RTGM, TARGET_RISK )

              
n = length( RTGMi );
switch n
    case 1
        COLORS = [ 'g' ];
    case 2
        COLORS = [ 'g' 'r' ];
    case 3
        COLORS = [ 'g' 'b' 'r' ];
    case 4
        COLORS = [ 'g' 'c' 'b' 'r' ];
    case 5
        COLORS = [ 'g' 'c' 'b' 'm' 'r' ];
    otherwise
        for i = 1:n, COLORS(i) = 'k';, end
        COLORS(1) = 'g';
        COLORS(end) = 'r';
end    

xmax = max( HazardCurve.SAs(end), max(RTGMi) ) / UHGM;
exponent = floor( abs(log10(xmax)) ) * sign( log10(xmax) );
XMAX = ceil( xmax / 10^exponent ) * 10^exponent;

xmin = min( HazardCurve.SAs(1), min(RTGMi) ) / UHGM;
exponent = ceil( abs(log10(xmin)) ) * sign( log10(xmin) );
XMIN = floor( xmin / 10^exponent ) * 10^exponent;

figure
set( gcf, 'Position', [ 9 49 944 948 ] )

% Hazard Curve
% ------------
subplot( 12, 1, 1:4 )
loglog( HazardCurve.SAs/UHGM, HazardCurve.AFEs/AFE4UHGM, 'k.-' )
xlim( [ XMIN XMAX ] )
hold on
loglog( HazardCurve.SAs(HazardCurve.Iextrap)/UHGM, ...
        HazardCurve.AFEs(HazardCurve.Iextrap)/AFE4UHGM, 'ko' )
loglog( xlim, ones(1,2), [ COLORS(1) ':' ], 'LineWidth', 2 )
for i = 1:n
    loglog( RTGMi(i)/UHGM*ones(1,2), ylim, [ COLORS(i) ':' ], 'LineWidth', 2 )
end
set( gca, 'XTickLabel', [] )
title( strvcat( 'RTGM Calculation Iterations', ' ' ), 'FontSize', 14 )
ylabel( strvcat( 'AFEs / AFE4UHGM', ' ' ), 'FontSize', 12 )

% Fragility PDF's
% ---------------
subplot( 12, 1, 5:6 )
for i = 1:n
    semilogx( FragilityCurves(i).SAs/UHGM, FragilityCurves(i).PDF, [ COLORS(i) '.-' ] )
    hold on
end
xlim( [ XMIN XMAX ] )
set( gca, 'XTickLabel', [] )
ylabel( strvcat( 'Fragility PDF', ' ' ), 'FontSize', 12 )

% Fragility CDF's
% ---------------
subplot( 12, 1, 7:8 )
for i = 1:n
    semilogx( FragilityCurves(i).SAs/UHGM, FragilityCurves(i).CDF, [ COLORS(i) '.-' ] )
    hold on
end
for i = 1:n
    semilogx( RTGMi(i)/UHGM*ones(1,2), ylim, [ COLORS(i) ':' ], 'LineWidth', 2 )
end
xlim( [ XMIN XMAX ] )
semilogx( xlim, FRAGILITY_AT_RTGM*ones(1,2), 'k:', 'LineWidth', 2 )
set( gca, 'XTickLabel', [] )
ylabel( strvcat( 'Fragility CDF', ' ' ), 'FontSize', 12 )

% Risk Integrands
% ---------------
subplot( 12, 1, 9:10 )
for i = 1:n
    RiskIntegrands(:,i) = FragilityCurves(i).PDF.*HazardCurve.AFEs;
    semilogx( FragilityCurves(i).SAs/UHGM, RiskIntegrands(:,i), [ COLORS(i) '.-' ] )
    hold on
end
xlim( [ XMIN XMAX ] )
set( gca, 'XTickLabel', [] )
ylabel( strvcat( 'Risk Integrand', ' ' ), 'FontSize', 12 )

% Risk Accumulations
% ------------------
subplot( 12, 1, 11:12 )
for i = 1:n
    RiskAccumulations(:,i) = cumtrapz( HazardCurve.SAs, RiskIntegrands(:,i) );
    semilogx( FragilityCurves(i).SAs/UHGM, RiskAccumulations(:,i) / TARGET_RISK, [ COLORS(i) '.-' ] )
    hold on
end
xlim( [ XMIN XMAX ] )
semilogx( xlim, ones(1,2), [ COLORS(end) ':' ], 'LineWidth', 2 )
xlabel( strvcat( ' ', 'SAs / UHGM' ), 'FontSize', 12 )
ylabel( strvcat( 'Risk / Target', ' ' ), 'FontSize', 12 )

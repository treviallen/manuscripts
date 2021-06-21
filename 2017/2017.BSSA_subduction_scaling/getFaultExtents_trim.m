% ***************************************************************************
% Program: getFaultExtents_1Seg.m
%
% Gets fault dimensions from the "Finite-Source Rupture Model Database"
% (http://www.seismo.ethz.ch/srcmod/Homepage.html) by Martin Mai and writes
% ShakeMap "*_fault.txt" file
%
% Author: T. Allen (20070730)
%
%% *************************************************************************

files = dir('matfiles/*.mat');
mlvect = [];
mwvect = [];
olvect = [];
owvect = [];
outtxt = [];


% Minimum percentage of cell slip of the maximum cell slip
minslippc = 0.10;

for i = 1:length(files)
%for i = 1:1
    geoX = [];
    geoXmin = [];
    geoY = [];
    geoXmin = [];
    geoZ = [];

    load(['matfiles/',files(i).name]);
    % get variable
    [pathstr,faultStrctStr,ext]=fileparts(files(i).name);
    %fs = eval(faultStrctStr);

    deg2rad = pi / 180;

    %% get cell size
    if fs.invSEGM == 1
        sfw = abs(fs.invDzDx(1)); % cell size
        sfl = abs(fs.invDzDx(2));

        %% get parameters of segment
        dep2top1 = fs.srcZ2top;
        wl1 = fs.srcDimWL;
        stk1 = fs.srcAStke;
        dip1 = fs.srcDipAn * deg2rad;

        slip = fs.slipSPL;

        maxSlip = max(slip(:));
        pcSlip = slip / maxSlip;

        % find cells gt 0.1
        [wcells lcells] = size(pcSlip);
        lvect = [];
        wvect = [];
        for j = 1:wcells
            ind = find(pcSlip(j,:) >= minslippc);
            ltmp = (max(ind) - min(ind)+1) * sfl; % lenght of rows
            if length(ind > 1)
                lvect = [lvect ltmp];
            end
            % test alt 
            geoX = [geoX fs.geoLON(j,max(ind)) fs.geoLON(j,min(ind))];
            geoY = [geoY fs.geoLAT(j,max(ind)) fs.geoLAT(j,min(ind))];
        end
        meanl = nanmean(lvect);
        meanl = prctile(lvect, 75);

        for k = 1:lcells
            ind = find(pcSlip(:,k) >= minslippc);
            wtmp = (max(ind) - min(ind)+1) * sfw; % column widths
            if length(ind > 1)
                wvect = [wvect wtmp];
            end
            geoZ = [geoY fs.geoZ(max(ind),k) fs.geoZ(min(ind),k)];
        end
        meanw = nanmean(wvect);
        meanw = prctile(wvect, 75);

        % get vectors for average change
        %mlvect = [mlvect meanl];
        %mwvect = [mwvect meanw];
        %olvect = [olvect fs.srcDimWL(2)];
        %owvect = [owvect fs.srcDimWL(1)];
        meanlvect = meanl;
        meanwvect = meanw;

        % get fault centroid
        fshape = size(fs.geoLAT);
        % do width
        if fshape(1)/2 == round(fshape(1)/2)
            if fshape(2)/2 == round(fshape(2)/2)
                centlat = mean(mean([fs.geoLAT(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
                centlon = mean(mean([fs.geoLON(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
                centdep = mean(mean([fs.geoZ(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
            else
                centlat = mean([fs.geoLAT(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
                centlon = mean([fs.geoLON(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
                centdep = mean([fs.geoZ(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
            end
        else
            if fshape(2)/2 == round(fshape(2)/2)
                centlat = mean([fs.geoLAT(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
                centlon = mean([fs.geoLON(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
                centdep = mean([fs.geoZ(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
            else
                centlat = mean([fs.geoLAT(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
                centlon = mean([fs.geoLON(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
                centdep = mean([fs.geoZ(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
            end
        end
        
        

        h = figure(i);
        
        %view(fs.geoViewAng, 30);
        set(gca,'Position',[0.13 0.11 0.775 0.815]);
        view(fs.srcAStke+90, 30);
        try
            plot3(centlon, centlat, str2num(centdep),'ro');
        catch err
            plot3(centlon, centlat, centdep,'ro');
        end

        %% from centroid, get coordinates
        halfhw = fs.srcDimWL(1)/2 * cos(fs.srcDipAn*deg2rad);
        halfv  = fs.srcDimWL(1)/2 * sin(fs.srcDipAn*deg2rad);

        % top coords
        [ctlat, ctlon] =  reckon(centlat, centlon, km2deg(halfhw), fs.srcAStke-90);
        hold on;
        plot3(ctlon, ctlat, centdep-halfv,'bo');
        tdep = centdep-halfv;

        [tclat1 tclon1] = reckon(ctlat, ctlon, km2deg(fs.srcDimWL(2)/2), fs.srcAStke);
        [tclat2 tclon2] = reckon(ctlat, ctlon, km2deg(fs.srcDimWL(2)/2), fs.srcAStke-180);

        % bottom coords
        halfhw = fs.srcDimWL(1)/2 * cos(fs.srcDipAn*deg2rad);
        halfv  = fs.srcDimWL(1)/2 * sin(fs.srcDipAn*deg2rad);

        [cblat, cblon] =  reckon(centlat, centlon, km2deg(halfhw), fs.srcAStke+90);
        hold on;
        plot3(cblon, cblat, centdep+halfv,'go');
        bdep = centdep+halfv;
        %set(gca,'ZGrid','on');

        [bclat1 bclon1]  = reckon(cblat, cblon, km2deg(fs.srcDimWL(2)/2), fs.srcAStke-180);
        [bclat2 bclon2]  = reckon(cblat, cblon, km2deg(fs.srcDimWL(2)/2), fs.srcAStke);

        oflat = [tclat1 tclat2 bclat1 bclat2 tclat1];
        oflon = [tclon1 tclon2 bclon1 bclon2 tclon1];
        ofdep = [tdep tdep bdep bdep tdep];
        hold on;
        plot3(oflon,oflat,ofdep,'k-','linewidth',2);

        % get average loc of pcSlip > 0.25
        slipind = find(pcSlip >= 0.25);
        sliplat = mean(fs.geoLAT(slipind));
        sliplon = mean(fs.geoLON(slipind));
        slipdep = mean(fs.geoZ(slipind));
        
        % TEST:  alt get centroid !!!!!!!!!!!!!!!!!!!! 
%         sliplon = mean(geoX);
%         sliplat = mean(geoY);
%         slipdep = mean(geoZ);
        

        %% from epicentre, get coordinates
        %% now get dimensions of trimmed fault
        halfhw = meanw/2 * cos(fs.srcDipAn*deg2rad);
        halfv  = meanw/2 * sin(fs.srcDipAn*deg2rad);

        % top coords
        [ctlat, ctlon] =  reckon(sliplat, sliplon, km2deg(halfhw), fs.srcAStke-90);
        hold on;
        %plot3(ctlon, ctlat, fs.evDPT-halfv,'bo');
        %plot3(cblon, cblat, slipdep-halfv,'mo');
        %tdep = centdep-halfv;
        tdep = slipdep-halfv;

        [tclat1 tclon1] = reckon(ctlat, ctlon, km2deg(meanl/2), fs.srcAStke);
        [tclat2 tclon2] = reckon(ctlat, ctlon, km2deg(meanl/2), fs.srcAStke-180);

        % bottom coords
        [cblat, cblon] =  reckon(sliplat, sliplon, km2deg(halfhw), fs.srcAStke+90);
        hold on;
        %plot3(cblon, cblat, fs.evDPT+halfv,'mo');
        %plot3(cblon, cblat, slipdep+halfv,'mo');
        %bdep = centdep+halfv;
        bdep = slipdep+halfv;


        [bclat1 bclon1]  = reckon(cblat, cblon, km2deg(meanl/2), fs.srcAStke-180);
        [bclat2 bclon2]  = reckon(cblat, cblon, km2deg(meanl/2), fs.srcAStke);

        mflat = [tclat1 tclat2 bclat1 bclat2 tclat1];
        mflon = [tclon1 tclon2 bclon1 bclon2 tclon1];
        mfdep = [tdep tdep bdep bdep tdep];
        hold on;
        plot3(mflon,mflat,mfdep,'r-','linewidth',2);
        %disp(mfdep);

        % catch negative depths and slide top edge of modified fault down-dip
        topind = [1 2 5];
        if tdep < 0
            vshift = abs(tdep) + 2.0; % add 2 km depth
            hshift = vshift / tan(fs.srcDipAn*deg2rad);
            [mflat(topind) mflon(topind)] = reckon(mflat(topind),mflon(topind), ...
                                            km2deg(hshift),fs.srcAStke+90);
            mfdep(topind) = mfdep(topind) + vshift;
            hold on;
            plot3(mflon,mflat,mfdep,'g--','linewidth',2);
        end
        
        % catch depths GT max depth and slide bottom edge of modified fault down-dip
        botind = [3 4];
        if bdep > max(max(fs.geoZ))
            vshift = bdep - max(max(fs.geoZ));
            hshift = vshift / tan(fs.srcDipAn*deg2rad);
            [mflat(botind) mflon(botind)] = reckon(mflat(botind),mflon(botind), ...
                                            km2deg(hshift),fs.srcAStke-90);
            mfdep(botind) = mfdep(botind) - vshift;
            hold on;
            plot3(mflon,mflat,mfdep,'g-','linewidth',2);
        end
        axis vis3d;
        %set(gca,'ZGrid','on');
        %colormap(hot)

        hold on;
        view(fs.srcAStke+90, 30);
        surf(fs.geoLON, fs.geoLAT, fs.geoZ, pcSlip);
        hold on;
        view(fs.srcAStke+90, 30);
        surf(fs.geoLON, fs.geoLAT, fs.geoZ, pcSlip);
        set(gca,'ZDir','reverse');
        grid on;

        zlabel('Depth (km)');
        xlabel('Longitude');
        ylabel('Latitude');
        view(fs.srcAStke+90, 30);
        title([fs.evDAT,' ',fs.event]);
        set(gca,'Position',[0.13 0.11 0.775 0.815]);
        refresh(h);
        
        % output file
        faultarray = [mflat' mflon' mfdep'];
        dlmwrite(['shake_fault/',fs.evTAG(2:9), ...
                 '_',fs.evTAG,'_fault.txt'], faultarray,'delimiter','\t')

    elseif fs.invSEGM > 9 && fs.invSEGM <= 10 % temp - change when can parse multiseg
        % loop through subfaults
        meanlvect = [];
        meanwvect = [];
        for ii=1:fs.invSEGM
            try
                sfw = eval(['fs.seg',num2str(ii),'invDzDx(1)']); % cell size
                sfl = eval(['fs.seg',num2str(ii),'invDzDx(2)']);
            catch err
                sfw = abs(fs.invDzDx(1)); % cell size
                sfl = abs(fs.invDzDx(2));
            end

            fstmp = fs;
            slip = eval(['fs.seg',num2str(ii),'SLIP']);

            maxSlip = max(slip(:));
            pcSlip = slip / maxSlip;

            % find cells gt 0.05
            [wcells lcells] = size(pcSlip);
            lvect = [];
            wvect = [];
            for j = 1:wcells
                ind = find(pcSlip(j,:) >= minslippc);
                ltmp = (max(ind) - min(ind)+1) * sfl;
                if length(ind > 1)
                    lvect = [lvect ltmp];
                end
            end
            %meanl = abs(nanmean(lvect));
            disp([sfl sfw])
            meanl = prctile(lvect, 75);

            for k = 1:lcells
                ind = find(pcSlip(:,k) >= minslippc);
                wtmp = (max(ind) - min(ind)+1) * sfw;
                if length(ind > 1)
                    wvect = [wvect wtmp];
                end
            end
            %meanw = abs(nanmean(wvect));
            meanw = prctile(wvect, 75);

            meanlvect = [meanlvect meanl];
            meanwvect = [meanwvect meanw];

            % get fault centroid
            geoLAT = eval(['fs.seg',num2str(ii),'geoLAT']);
            geoLON = eval(['fs.seg',num2str(ii),'geoLON']);
            geoZ   = eval(['fs.seg',num2str(ii),'geoZ']);
            fshape = size(geoLAT);

            % do width
            if fshape(1)/2 == round(fshape(1)/2)
                if fshape(2)/2 == round(fshape(2)/2)
                    centlat = mean(mean([geoLAT(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
                    centlon = mean(mean([geoLON(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
                    centdep = mean(mean([geoZ(fshape(1)/2:fshape(1)/2+1, fshape(2)/2:fshape(2)/2+1)]));
                else
                    centlat = mean([geoLAT(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
                    centlon = mean([geoLON(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
                    centdep = mean([geoZ(fshape(1)/2:fshape(1)/2+1, ceil(fshape(2)/2))]);
                end
            else
                if fshape(2)/2 == round(fshape(2)/2)
                    centlat = mean([geoLAT(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
                    centlon = mean([geoLON(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
                    centdep = mean([geoZ(ceil(fshape(1)/2), fshape(2)/2:fshape(2)/2+1)]);
                else
                    centlat = mean([geoLAT(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
                    centlon = mean([geoLON(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
                    centdep = mean([geoZ(ceil(fshape(1)/2), ceil(fshape(2)/2))]);
                end
            end

            h = figure(i);
            set(gca,'Position',[0.13 0.11 0.775 0.815]);

            try
                plot3(centlon, centlat, str2num(centdep),'ro');
            catch err
                plot3(centlon, centlat, centdep,'ro');
            end

            %% from centroid, get coordinates
            segDipAn = eval(['fs.seg',num2str(ii),'DipAn']);
            segStkAn = eval(['fs.seg',num2str(ii),'AStke']);
            segDimWL = eval(['fs.seg',num2str(ii),'DimWL']);
            halfhw = segDimWL(1)/2 * cos(segDipAn*deg2rad);
            halfv  = segDimWL(1)/2 * sin(segDipAn*deg2rad);

            % top coords
            [ctlat, ctlon] =  reckon(centlat, centlon, km2deg(halfhw), segStkAn-90);
            hold on;
            plot3(ctlon, ctlat, centdep-halfv,'bo');
            tdep = centdep-halfv;

            [tclat1 tclon1] = reckon(ctlat, ctlon, km2deg(segDimWL(2)/2), segStkAn);
            [tclat2 tclon2] = reckon(ctlat, ctlon, km2deg(segDimWL(2)/2), segStkAn-180);

            % bottom coords
            halfhw = segDimWL(1)/2 * cos(segDipAn*deg2rad);
            halfv  = segDimWL(1)/2 * sin(segDipAn*deg2rad);

            [cblat, cblon] =  reckon(centlat, centlon, km2deg(halfhw), segStkAn+90);
            hold on;
            plot3(cblon, cblat, centdep+halfv,'go');
            bdep = centdep+halfv;
            %set(gca,'ZGrid','on');

            [bclat1 bclon1]  = reckon(cblat, cblon, km2deg(segDimWL(2)/2), segStkAn-180);
            [bclat2 bclon2]  = reckon(cblat, cblon, km2deg(segDimWL(2)/2), segStkAn);

            oflat = [tclat1 tclat2 bclat1 bclat2 tclat1];
            oflon = [tclon1 tclon2 bclon1 bclon2 tclon1];
            ofdep = [tdep tdep bdep bdep tdep];
            hold on;
            plot3(oflon,oflat,ofdep,'linewidth',2);

            % get average loc of pcSlip > 0.5
            slipind = find(pcSlip >= 0.25);
            sliplat = mean(geoLAT(slipind));
            sliplon = mean(geoLON(slipind));
            slipdep = mean(geoZ(slipind));

            %% from epicentre, get coordinates
            %% now get dipensions of trimmed fault
            halfhw = meanw/2 * cos(segDipAn*deg2rad);
            halfv  = meanw/2 * sin(segDipAn*deg2rad);

            % top coords
            [ctlat, ctlon] =  reckon(sliplat, sliplon, km2deg(halfhw), segStkAn-90);
            hold on;
            %plot3(ctlon, ctlat, fs.evDPT-halfv,'bo');
            %plot3(cblon, cblat, slipdep-halfv,'mo');
            %tdep = centdep-halfv;
            tdep = slipdep-halfv;

            [tclat1 tclon1] = reckon(ctlat, ctlon, km2deg(meanl/2), segStkAn);
            [tclat2 tclon2] = reckon(ctlat, ctlon, km2deg(meanl/2), segStkAn-180);

            % bottom coords
            [cblat, cblon] =  reckon(sliplat, sliplon, km2deg(halfhw), segStkAn+90);
            hold on;
            %plot3(cblon, cblat, fs.evDPT+halfv,'mo');
            %plot3(cblon, cblat, slipdep+halfv,'mo');
            %bdep = centdep+halfv;
            bdep = slipdep+halfv;


            [bclat1 bclon1]  = reckon(cblat, cblon, km2deg(meanl/2), segStkAn-180);
            [bclat2 bclon2]  = reckon(cblat, cblon, km2deg(meanl/2), segStkAn);

            mflat = [tclat1 tclat2 bclat1 bclat2 tclat1];
            mflon = [tclon1 tclon2 bclon1 bclon2 tclon1];
            mfdep = [tdep tdep bdep bdep tdep];
            hold on;
            plot3(mflon,mflat,mfdep,'r-','linewidth',2);
            %disp(mfdep);

            % catch negative depths and slide modified fault down-dip
            if tdep < 0
                vshift = abs(tdep) + 2;
                hshift = vshift / tan(segDipAn*deg2rad);
                [mflat(topind) mflon(topind)] = reckon(mflat(topind),mflon(topind), ...
                                                km2deg(hshift),fs.srcAStke+90);
                mfdep(topind) = mfdep(topind) + vshift;
                hold on;
                plot3(mflon,mflat,mfdep,'g-','linewidth',2);
            end
            axis vis3d;
            %set(gca,'ZGrid','on');
            %colormap(hot)

            hold on;
            view(segStkAn+90, 30);
            surf(geoLON, geoLAT, geoZ, pcSlip);
            set(gca,'ZDir','reverse');
            grid on;

            zlabel('Depth (km)');
            xlabel('Longitude');
            ylabel('Latitude');
            view(segStkAn+90, 30);
            title([fs.evDAT,' ',fs.event]);
            set(gca,'Position',[0.13 0.11 0.775 0.815]);
            refresh(h);
       end


    end
    mlvect = [mlvect sum(meanlvect)];
    mwvect = [mwvect mean(meanwvect)];
    olvect = [olvect fs.srcDimWL(2)];
    owvect = [owvect fs.srcDimWL(1)];

    % make txt
%     outtxt = [outtxt datestr(datenum(fs.evDAT,'yyyy/mm/dd'),'yyyymmdd'),',', ...
%               num2str(fs.srcDimWL(2)),',',num2str(fs.srcDimWL(1)),',', ...
%               num2str(sum(meanlvect)),',',num2str(mean(meanwvect)) char(10)];
          
    outtxt = [outtxt fs.evTAG(2:13),',', num2str(fs.srcMwMoS(1)),',', ...
              num2str(fs.evLON),',',num2str(fs.evLAT),',',num2str(fs.evDPT),',', ...
              num2str(fs.srcDimWL(2)),',',num2str(fs.srcDimWL(1)),',', ...
              num2str(sum(meanlvect)),',',num2str(mean(meanwvect)) char(10)];     

    % save figure
    %saveas(gcf,[datestr(datenum(fs.evDAT,'mm/dd/yyyy'),'yyyymmdd'),'.jpg'],'jpg');
    saveas(gcf,['shake_fault/',fs.evTAG(2:9),'.jpg'],'jpg');
end

% get percentage differences
pcldiff = mlvect ./ olvect;
pcwdiff = mwvect ./ owvect;

header = 'DATE,MW,LON,LAT,DEP,OLEN,OWID,MLEN,MWID';
outtxt = [header,char(10),outtxt];
dlmwrite(trimmed_faults.csv',outtxt,'delimiter','');
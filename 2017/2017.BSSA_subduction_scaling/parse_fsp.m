

fspfiles = dir('FSPs/c*.fsp');
%fspfile = 'FSPs/c201103110546.fsp';

for i = 1:length(fspfiles)

%for i = 1:1
    clear -v fs

    % load fsp file
    fspfile = ['FSPs/',fspfiles(i).name];
    fid = fopen(fspfile);
    fgetl(fid); fgetl(fid);

    % get event name and date
    line = fgetl(fid);
    %event = strsplit(line,'  ');
    %evdate = strsplit(event{2},' ');
    event = strsplit('  ',line);
    evdate = strsplit(' ',event{2});
    fs.evDAT = evdate{1};
    %evdate = strsplit(evdate{1}, '/')

    event = regexprep(event(1),'% Event : ','');
    fs.event = event{1};

    % get event tag
    line = fgetl(fid);
    fs.evTAG = regexprep(line, '% EventTAG: ','');

    % get event location
    fgetl(fid);
    line = fgetl(fid);
    %loc = strsplit(line,' ');
    loc = strsplit(' ',line);
    fs.evLAT = str2double(loc{7});
    fs.evLON = str2double(loc{11});
    fs.evDPT = str2double(loc{15});

    % get size
    line = fgetl(fid);
    %siz = strsplit(line,' ');
    siz = strsplit(' ',line);
    fs.srcDimWL = [str2double(siz{11}) str2double(siz{6})];
    fs.srcMwMoS = [str2double(siz{16}) str2double(siz{20})];

    % get strike, dip, rake
    line = fgetl(fid);
    %sdr = strsplit(line,' ');
    sdr = strsplit(' ',line);
    fs.srcARake = str2double(sdr{14});
    fs.srcAStke = str2double(sdr{6});
    fs.srcDipAn = str2double(sdr{10});

    % get inversion constraints
    fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
    line = fgetl(fid);
    %inv = strsplit(line,' ');
    inv = strsplit(' ',line);
    fs.invNzNx =  [str2double(inv{10}) str2double(inv{6})];
    line = fgetl(fid);
    %inv = strsplit(line,' ');
    inv = strsplit(' ',line);
    fs.invDzDx =  [str2double(inv{11}) str2double(inv{6})];

    % get num segs
    line = fgetl(fid);
    %seg = strsplit(line,' ');
    seg = strsplit(' ',line);
    fs.invSEGM = str2double(seg{10});

    % now have all header info, get fault data

    if fs.invSEGM == 1

        % get number of headerlines
        j = 0;
        while strcmp(line,'% LAT LON X==EW Y==NS Z SLIP RAKE TRUP RISE') == 0
            line = fgetl(fid);
            j = j+1;
        end
        [lat lon ew ns z slip rake trup rise] = textread(fspfile, '%f%f%f%f%f%f%f%f%f', ...
                                                         'headerlines',16+j);


        % make rupture matrices
        if z(1) == z(2)
            fs.geoLAT = zeros(fs.invNzNx(1),fs.invNzNx(2));
            fs.geoLON = zeros(fs.invNzNx(1),fs.invNzNx(2));
            fs.geoX = zeros(fs.invNzNx(1),fs.invNzNx(2));
            fs.geoY = zeros(fs.invNzNx(1),fs.invNzNx(2));
            fs.geoZ = zeros(fs.invNzNx(1),fs.invNzNx(2));
            fs.slipSPL = zeros(fs.invNzNx(1),fs.invNzNx(2));
            for x = 1:fs.invNzNx(2)
                fs.geoLAT(:,x) = lat(x:fs.invNzNx(2):end);
                fs.geoLON(:,x) =  lon(x:fs.invNzNx(2):end);
                fs.geoX(:,x) =    ew(x:fs.invNzNx(2):end);
                fs.geoY(:,x) =    ns(x:fs.invNzNx(2):end);
                fs.geoZ(:,x) =    z(x:fs.invNzNx(2):end);
                fs.slipSPL(:,x) = slip(x:fs.invNzNx(2):end);
            end
        else
            fs.geoLAT = [];
            fs.geoLON = [];
            fs.geoX = [];
            fs.geoY = [];
            fs.geoZ = [];
            fs.slipSPL = [];
            for x = 1:fs.invNzNx(2)
                nstart = (x-1)*fs.invNzNx(1)+1;
                fs.geoLAT =  [fs.geoLAT lat(nstart:nstart+fs.invNzNx(1)-1)];
                fs.geoLON =  [fs.geoLON lon(nstart:nstart+fs.invNzNx(1)-1)];
                fs.geoX =    [fs.geoX ew(nstart:nstart+fs.invNzNx(1)-1)];
                fs.geoY =    [fs.geoY ns(nstart:nstart+fs.invNzNx(1)-1)];
                fs.geoZ =    [fs.geoZ z(nstart:nstart+fs.invNzNx(1)-1)];
                fs.slipSPL = [fs.slipSPL slip(nstart:nstart+fs.invNzNx(1)-1)];
            end
        end
        fs.srcZ2top = min(z);

    elseif fs.invSEGM > 1

        % get number of headerlines
        j = 0;
        while strcmp(line,'%--------------------------- MULTISEGMENT MODEL ---------------------------------------------------') == 0
            line = fgetl(fid);
            j = j+1;
        end
        line = fgetl(fid);
        
        % now loop through subfaults
        for k = 1:fs.invSEGM
            line = fgetl(fid);
            sdr = strsplit(' ',line);
            %sdr = strsplit(line,' ');
            eval(['fs.seg',num2str(k),'AStke = str2double(sdr{7})']);
            eval(['fs.seg',num2str(k),'DipAn = str2double(sdr{11})']);

            line = fgetl(fid);
            wl = strsplit(' ',line);
            %wl = strsplit(line,' ');
            eval(['fs.seg',num2str(k),'DimWL = [str2double(wl{8}) str2double(wl{4})]']);

            line = fgetl(fid);
            z2t = strsplit(' ',line);
            %z2t = strsplit(line,' ');
            eval(['fs.seg',num2str(k),'Z2top = str2double(z2t{7})']);

            % read lines
            fgetl(fid); fgetl(fid); fgetl(fid); line = fgetl(fid);
            sf = strsplit(' ',line);
            %sf = strsplit(line,' ');
            nsubf = str2double(sf{4});

            % get fault data
            [lat lon ew ns z slip rake trup rise] = textread(fspfile, '%f%f%f%f%f%f%f%f%f', nsubf, ...
                                                             'headerlines',26+j);
            j = j + nsubf + 11;

            % get segment NzNx
              eval(['fs.seg',num2str(k),'NzNx = [round(fs.seg',num2str(k),'DimWL(1) / fs.invDzDx(1)) round(fs.seg', ...
                    num2str(k),'DimWL(2) / fs.invDzDx(2))]']);

            % make rupture matrices
            if z(1) == z(2)
                eval(['fs.seg',num2str(k),'geoLAT  = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                eval(['fs.seg',num2str(k),'geoLON  = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                eval(['fs.seg',num2str(k),'geoX    = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                eval(['fs.seg',num2str(k),'geoY    = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                eval(['fs.seg',num2str(k),'geoZ    = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                eval(['fs.seg',num2str(k),'slipSPL = zeros(fs.seg',num2str(k),'NzNx(1),fs.seg',num2str(k),'NzNx(2))']);
                for x = 1:eval(['fs.seg',num2str(k),'NzNx(2)'])
                    eval(['fs.seg',num2str(k),'geoLAT(:,x)  = lat(x:fs.seg',num2str(k),'NzNx(2):end)']);
                    eval(['fs.seg',num2str(k),'geoLON(:,x)  = lon(x:fs.seg',num2str(k),'NzNx(2):end)']);
                    eval(['fs.seg',num2str(k),'geoX(:,x)    =  ew(x:fs.seg',num2str(k),'NzNx(2):end)']);
                    eval(['fs.seg',num2str(k),'geoY(:,x)    =  ns(x:fs.seg',num2str(k),'NzNx(2):end)']);
                    eval(['fs.seg',num2str(k),'geoZ(:,x)    =   z(x:fs.seg',num2str(k),'NzNx(2):end)']);
                    eval(['fs.seg',num2str(k),'slipSPL(:,x) =slip(x:fs.seg',num2str(k),'NzNx(2):end)']);
                end
            else
                fs.geoLAT = [];
                fs.geoLON = [];
                fs.geoX = [];
                fs.geoY = [];
                fs.geoZ = [];
                fs.slipSPL = [];
                for x = 1:fs.invNzNx(2)
                    nstart = (x-1)*fs.invNzNx(1)+1;
                    fs.geoLAT =  [fs.geoLAT lat(nstart:nstart+fs.invNzNx(1)-1)];
                    fs.geoLON =  [fs.geoLON lon(nstart:nstart+fs.invNzNx(1)-1)];
                    fs.geoX =    [fs.geoX ew(nstart:nstart+fs.invNzNx(1)-1)];
                    fs.geoY =    [fs.geoY ns(nstart:nstart+fs.invNzNx(1)-1)];
                    fs.geoZ =    [fs.geoZ z(nstart:nstart+fs.invNzNx(1)-1)];
                    fs.slipSPL = [fs.slipSPL slip(nstart:nstart+fs.invNzNx(1)-1)];
                end
            end
            % read to next subfault
            for n = 1:nsubf+4
                fgetl(fid);
            end
            fs.srcZ2top = min(z);
        end
    end
    % save to mat file
    filename = ['matfiles/',fs.evTAG,'.mat'];
    save(filename,'fs');

    %disp(fs)
    fclose(fid);
end
- gmt5 surface hazard_map-mean_PGA_0.02.csv -Gnsha18_0.02_interp.0.05.grd -R110/156/-46/-9 -I0.05
- gmt5 grdmath nsha18_0.02_interp.0.05.grd 0.002 MAX = nsha18_0.02_interp.0.05.grd
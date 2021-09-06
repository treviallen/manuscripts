* gmt5 grd2xyz gshap_au.grd > gshap.xyz
* gmt5 surface gshap.xyz -Ggshap_interp.0.05.grd -R110/156/-46/-9 -I0.05
* gmt5 grdmath gshap_interp.0.05.grd 9.80665 DIV = gshap_interp.0.05.grd
The SMSIM code has been modified to run in commandline and not pause for input control file.  Instead, the source code has been modified such that the control file "allen07_a_ts_drvr.ctl" will be called by default.  This file is produced dynamically by the code "simulate_acceleration_data.py", i.e.:

	calling "./a_ts_drvr_mod" is the same as calling "./a_ts_drvr allen07_a_ts_drvr.ctl"
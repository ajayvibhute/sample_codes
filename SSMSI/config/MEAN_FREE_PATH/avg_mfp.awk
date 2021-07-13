
####################################################################
#            Mission: ASTROSAT
#            Project: SSM
#
#          File-Name: avg_mfp.awk ($Revision: 1.2 $)
#%%             Path:
#%% /data1/ravibt/ASTROSAT/SSM/etc/MEAN_FREE_PATH
#         Created on: 2008 Jul 22 01:51:19 PM (IST: UTC+0530)
#   Last Modified on: 2008 Jul 25 10:26:14 AM (IST: UTC+0530)
#             Author: Ravishankar B.T. (ravibt@isac.gov.in)
#        Supervisors: Dipankar Bhattacharya (dipankar@iucaa.ernet.in)
#                     Seetha S (seetha@isac.gov.in)
####################################################################



####################################################################
#
# $Id: avg_mfp.awk,v 1.2 2008/07/25 05:53:54 ravibt Exp $
# $Log: avg_mfp.awk,v $
# Revision 1.2  2008/07/25 05:53:54  ravibt
# Will not do a "flat" average, but spectrum slope weighted average.
#
# Revision 1.1  2008/07/23 06:27:11  ravibt
# Initial revision
#
####################################################################


####################################################################
#
# Usage: awk -f avg_mfp.awk mfp.dat > avg_mfp.dat
#
# The first two columns of input file ("mfp.dat") are expected to be
# (1) Energy, keV and (2) Mean-Free-path (in any units - output will
# also be in these units).
#
# Initial Revision: Assuming a "flat" average of Mean-free-paths.
# Revision 1.1: Spectrum slope weighted average of Mean-free-paths.
#               (for now assuming Crab spectra - N(E) ~ E{-2})
#
####################################################################


{
	if($1 != "#"){

		if($1 >= 2.7 && $1 < 4.8){
			band1_wt_mfp += $2/($1*$1)
			band1_wt += 1/($1*$1)
			nband1++;
		}

		if($1 >= 4.8 && $1 < 6.5){
			band2_wt_mfp += $2/($1*$1)
			band2_wt += 1/($1*$1)
			nband2++;
		}

		if($1 >= 6.5 && $1 <= 10){
			band3_wt_mfp += $2/($1*$1)
			band3_wt += 1/($1*$1)
			nband3++;
		}

	}
}
END{
	print band1_wt_mfp/band1_wt "\t" band2_wt_mfp/band2_wt "\t" band3_wt_mfp/band3_wt
}

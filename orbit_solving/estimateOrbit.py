import numpy as np
import matplotlib.pyplot as plt
import argparse
#Just plot baricentric period vs epoch.
#sys.argv[1] should be in the format of {pulsar_name}.per.
#sys.argv[2] should be the range of ephemeris you want to plot, as in start:end, staring with 0.
#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 08/09/2021

def fit_ellipse(periods,derivatives):

	print("Assuming circular orbit and guessing...")
	print("")

	coeff=np.polynomial.polynomial.polyfit(history[1],history[2]**2,2)

	period=-coeff[1]/(2*coeff[2])
	p_orb=2*np.pi*299792458/(period*np.sqrt(-coeff[2]))
	axis=p_orb*np.sqrt(period**2-coeff[0]/coeff[2])/(2*np.pi*period)
	acc_amp=np.square(2*np.pi/p_orb)*axis*299792458
	p_amp=2*np.pi*period*axis/p_orb

	print("Spin period= {} ms".format(period*1000))
	print("Orbital period= {} days".format(p_orb/24/3600))
	print("Projected axis= {} ls".format(axis))
	phase=np.arctan(-history[2]/(acc_amp*(history[1]-period)*history[1]))
	per_epoch=history[0]-phase*p_orb/(2*np.pi*24*3600)
	print("First epoch of periastron= {} MJD".format(per_epoch[0]))
	print("")

	return (period,p_orb,axis,acc_amp,p_amp,coeff)

def fit_folding(epochs,periods,folding_range="none"):

	if folding_range != "none":
		trial_porb=float(folding_range.split(":")[0])
		interval=float(folding_range.split(":")[1])
	else:
		interval=epochs[-1]-epochs[0]
		trial_porb=0.041
	smallest_interval=np.min(epochs[1:]-epochs[:-1])
	epochs=epochs-epochs[0]

	print("Searching for best orbital period in between {} and {} days...".format(trial_porb,interval))
	print("")

	roughness=[]
	trial_porbs=[]
	i=0

	while trial_porb <= interval:
		if i/10000==i//10000:
			print("Current period: {} days. Current step: {} days".format(trial_porb,[0.1*smallest_interval*trial_porb/interval,1e-2*(trial_porb**2)/(2*np.pi*interval)]))
		folded_epochs=epochs-(epochs//trial_porb)*trial_porb
		folded_periods=periods[np.argsort(folded_epochs)]
		roughness.append(np.sum((folded_periods-np.roll(folded_periods,-1))**2))
		trial_porbs.append(trial_porb)
		trial_porb=trial_porb+min(0.1*smallest_interval*trial_porb/interval,1e-2*(trial_porb**2)/(2*np.pi*interval))
#		trial_porb=trial_porb+1e-5*(trial_porb**2)/(2*np.pi*interval)
		i=i+1

	roughness=np.array(roughness)
	trial_porbs=np.array(trial_porbs)
	print("Total number of trials: {}".format(np.size(roughness)))

	best_porb=trial_porbs[np.argmin(roughness)]

	print("Best orbital period= {} days".format(best_porb))
	print("")

	return (best_porb,trial_porbs,roughness)

parser=argparse.ArgumentParser()
parser.add_argument("data",help="File with columns of MJD, period (and derivatives, derivative uncertainty).")
#parser.add_argument("-e","--ephemeris",help="tempo2 file with an ephemeris.")
parser.add_argument("-p","--period",help="Units of period. Default: 'ms'.",choices=["s","ms","s-1"])
parser.add_argument("-d","--derivative",help="Units of derivative. Default: 's/s'.",choices=["s/s","s-2","m/s2"])
parser.add_argument("-r","--range",help="Range of data lines to be read 'min:max'. Default: '0:inf'.")
parser.add_argument("-m","--method",help="Method to fit the data with",choices=["ellipse","folding"])
parser.add_argument("--fold_range",help="Custom orbital period range in serch with -m='folding' in days. 'min:max' in days")
parser.add_argument("-M","--methods",help="Print details about methods.",action="store_true")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

if args.methods:
	print("Currently, there are 2 models t fit the data with:")
	print("")
	print(" - Ellipse:")
	print("   It reads the 2nd (Period), 3rd (Derivative) and 4th (Derivative uncertainty) columns of data and apply Freire, Kramer & Lyne (2000). It fits a parabola on the absolute values of acceleration in function of baricentric period. The coefficients tell you about the orbit. Very good if you have good period and derivative values and the orbit is not excentric.")
	print("")
	print(" - Folding:")
	print("   It reads the 1st (MJD) and 2nd (Period) columns of data and apply Bhattacharyya & Nityananda (2008). It folds your data points at consecutive trial periods in between 1 hour and the duration of the data set, estimating a roughness value that should be minimal at the true period. It can take a long time to run, but is is very useful if the orbit is eccentric or you don't have good derivative measureements, but requires many data points on a cadence higher than the orbital period to work.")
	print("   If '--fold_range' is specified, then the orbital period search range is restricted in days to the custom values.")
	exit()

if args.period:
	p_units=args.period
else:
	p_units="ms"
if args.verbose==True:
	print("Period units are "+p_units+".")
	print("")

if args.derivative:
	pdot_units=args.derivative
else:
	pdot_units="s/s"
if args.verbose==True:
	print("Derivative units are "+pdot_units+".")
	print("")

if args.range:
	start=int(args.args.range.split(":")[0])
	end=int(args.range.split(":")[1])
	history=np.loadtxt(data).T[:,start:end+1]
	if args.verbose==True:
		print("Reading ephemeris from {} to {} from file {}".format(start,end,args.data))
		print("")
else:
	history=np.loadtxt(args.data).T
	if args.verbose==True:
		print("Reading all ephemeris from file {}".format(args.data))
		print("")

if args.method:
	fit=True
	model=args.method
	if args.fold_range:
		folding_range=args.fold_range
	if args.verbose:
		print("The data will be fitted with the followng model: "+args.method)
		if model=="folding" and args.fold_range:
			print("The folding range will be of "+folding_range+" days.")
		print("Write '-M' to know more about the models.")
		print("")

else:
	fit=False
	model="none"
	if args.verbose:
		print("The data will just be show in a period - epoch diagram")
		print("Give '-m' to estimate orbit with a model.")
		print("Write '-M' to know more about the models.")
		print("")


if p_units=="s-1":
	history[1]=1/history[1]

if p_units=="ms":
	history[1]=history[1]/1000

if pdot_units=="s/s" and model=="ellipse":
	history[2]=299792458*history[2]/history[1]
	history[3]=299792458*history[3]/history[1]

if pdot_units=="s-2" and model=="ellipse":
	history[2]=299792458*(-history[2]*(history[1]**2))/history[1]
	history[3]=299792458*(-history[3]*(history[1]**2))/history[1]

if fit==True:
	if model == "ellipse":
		(period,p_orb,axis,acc_amp,p_amp,coeff)=fit_ellipse(history[1],history[2])
	elif model == "folding":
		if 'folding_range' in locals():
			(p_orb,trial_porbs,roughness)=fit_folding(history[0],history[1],folding_range=folding_range)
		else:
			(p_orb,trial_porbs,roughness)=fit_folding(history[0],history[1])
	else:
		sys.exit("Please state a valid orbit estimation algorithm.")

if fit == False:

	plt.plot(history[0]-int(history[0,0]),history[1]*1000,"o")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("epoch - {} (MJD)".format(int(history[0,0])))
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.show()

#	plt.errorbar(1000*history[1],history[2]**2,yerr=2*history[2]*history[3],fmt="o")
#	plt.ylabel("$acc^2$ (m/s$^2$)$^2$")
#	plt.xlabel("$P_{bary}$ (ms)")
#	plt.title(sys.argv[-1].split(".")[0])
#	plt.tight_layout()
#	plt.show()

#	plt.errorbar(1000*history[1],history[2],yerr=history[3],fmt="o")
#	plt.ylabel("$acc$ (m/s$^2$)")
#	plt.xlabel("$P_{bary}$ (ms)")
#	plt.title(sys.argv[-1].split(".")[0])
#	plt.tight_layout()
#	plt.legend()
#	plt.show()

elif fit == True and model == "ellipse":
	
	plt.errorbar(1000*history[1],history[2]**2,yerr=2*history[2]*history[3],fmt="o")
	periods=np.linspace(period-p_amp,period+p_amp,1000)
	plt.plot(1000*periods,coeff[0]+coeff[1]*periods+coeff[2]*periods**2,"c-")
	plt.ylabel("$acc^2$ (m/s$^2$)$^2$")
	plt.xlabel("$P_{bary}$ (ms)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.show()

	plt.errorbar(1000*history[1],history[2],yerr=history[3],fmt="o")
	phases=np.linspace(0,2*np.pi,1000)
	periods=period+p_amp*np.cos(phases)
	accs=-acc_amp*np.sin(phases)
	plt.plot(1000*periods,accs,"c-",label="$P_{orb}=$ "+str(round(p_orb/(24*3600),2))+" d, $x=$ "+str(round(axis,2))+" ls")
	plt.ylabel("$acc$ (m/s$^2$)")
	plt.xlabel("$P_{bary}$ (ms)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()

	epochs=history[0,:]-history[0,0]
	folded_epochs=epochs-(epochs//(p_orb/(24*3600)))*(p_orb/(24*3600))
	#To do: Add lines to also plot a sinusoid with the orbit, perhaps.
	plt.plot(folded_epochs,history[1]*1000,"o",label="$P_{orb}=$ "+str(round(p_orb/(24*3600),2))+" d")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("orbital phase (MJD)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()

elif fit == True and model == "folding":
	
	plt.plot(trial_porbs,roughness/np.max(roughness))
	plt.ylabel("Normalized $Roughness$")
	plt.xlabel("Trial orbital period (days)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.show()

	epochs=history[0,:]-history[0,0]
	folded_epochs=epochs-(epochs//p_orb)*p_orb
	plt.plot(folded_epochs,history[1]*1000,"o",label="$P_{orb}=$ "+str(round(p_orb,5))+" d")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("orbital phase (MJD)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()


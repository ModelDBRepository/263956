#======================================================
# Simulation parameters #
#======================================================

IC=1
NumSim=1
ssnn=121

N=1025
T=8000

m=0 #hierarchical level

type_inh=2	# 1 = FS ; 2 = LTS; 3 = RS
type_ex1=1 # 1 = RS ; 2 = IB ; 3 = CH
type_ex2=3	# 1 = RS ; 2 = IB ; 3 = CH

perEx2=20

#======================================================
sgin=100
sgex=15

Noise_value=10000 #=1E-5
noiseAmp=`echo "scale=10;$Noise_value / 1000000000" | bc`

echo $noiseAmp
dtStepsTraj=0
NameFileIC="raster.dat"

perc_exc1=`echo "scale=1;(100- $perEx2)" | bc`
perc_exc2=`echo "scale=1;$perEx2" | bc`
gin=`echo "scale=3; $sgin / 100" | bc`
gex=`echo "scale=3; $sgex / 100" | bc`

#======================================================
# Parameters to create stimulation #
#======================================================

in_stim=0 
InpAmp=0
tStim=0

TimeFreeFall=0
dtTraj=0

changeNneurons=0
tStimAfterChange=0
InpAmpAfterChange=0

#======================================================
# Run #
#======================================================

time ./code.out $m $perc_exc1 $perc_exc2 $tStim $InpAmp $in_stim $gin $gex $ssnn $NameFileIC $type_inh $type_ex1 $type_ex2 $N $T $NumSim $IC $changeNneurons $tStimAfterChange $InpAmpAfterChange $TimeFreeFall $dtTraj $dtStepsTraj $noiseAmp &
		


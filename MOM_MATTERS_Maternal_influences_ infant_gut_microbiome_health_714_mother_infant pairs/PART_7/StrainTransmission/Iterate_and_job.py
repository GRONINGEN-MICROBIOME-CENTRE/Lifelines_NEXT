from pathlib import Path
import subprocess
import time
def check_job_count():
	output = subprocess.check_output(['squeue', '-u', 'umcg-sandreusanchez']).decode('utf-8')
	job_count = len(output.splitlines()) - 1
	return(job_count)

def Script1_transmission():
	script = "1_Assess_transmission.R"

	for File in Path("NEXT_RAxML_distmats").glob("*.txt"):
		SGB = File.stem.replace("_DistMat", "")
		if Path("Distance_tables/"+SGB+".rds").exists(): continue
		Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))	
		Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.StrainTransm
#SBATCH --output=logs/{SGB}.StrainTransm.o
#SBATCH --error=logs/{SGB}.StrainTransm.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:20:00
#SBATCH --mem=20GB

cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
		Script_n = "scripts/{SGB}_StrainTransm.sh".format(SGB = SGB)
		with open(Script_n,'w') as F: F.write(Job_setup + Command)

		job_count = check_job_count()
		while job_count >= 100:
			time.sleep(300)
			job_count = check_job_count()
		subprocess.call("sbatch "+Script_n, shell=True)
		#print(Script_n)
def Script2_LinearModel_transmission():
	SGBs_run = []
	script = "2_MotherInfantTranmission_Association.R"
	with open("SGB_babies.tsv") as Info:
		for line in Info:
			l = line.rstrip().split()
			if l[0] == "SGB": continue
			if int(l[1]) < 20: continue
			SGBs_run.append(l[0])
	SGBs_run = ["t__SGB17248", "t__SGB1861"]
	for SGB in SGBs_run:
		Out = "Results/{S}_summaryLMERTEST_COV_infanttime+mothertime+nextid.tsv".format(S = SGB)
		#if Path(Out).exists(): continue
		#call second script
		File = "Distance_tables/{S}.rds".format(S=SGB)
		Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))
		Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.LM
#SBATCH --output=logs/{SGB}.LM.o
#SBATCH --error=logs/{SGB}.LM.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=20GB
cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
		Script_n = "scripts/{SGB}_StrainTransmLM.sh".format(SGB = SGB)
		with open(Script_n,'w') as F: F.write(Job_setup + Command)
		subprocess.call("sbatch "+Script_n, shell=True)
		#print(Script_n)

def Script3_mantelTest():
	script = "4_MantelTest_Geography_vs_Distance.R"
	for File in Path("NEXT_RAxML_distmats").glob("*.txt"):
		SGB = File.stem.replace("_DistMat", "")
		if Path("Mantel_distance/",SGB,".tsv").exists(): continue
		Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))
		Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.mantel
#SBATCH --output=logs/{SGB}.mantel.o
#SBATCH --error=logs/{SGB}.mantel.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:50:00

cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
		Script_n = "scripts/{SGB}_mantel.sh".format(SGB = SGB)
		with open(Script_n,'w') as F: F.write(Job_setup + Command)
		job_count = check_job_count()
		while job_count >= 100:
			time.sleep(60)
			job_count = check_job_count()
		subprocess.call("sbatch "+Script_n, shell=True)
		#print(Script_n)

def Secript_motherConsistency():
	script = "6_MaternalStrainStability.R"
	for File in Path("Distance_tables").glob("*.rds"):
		SGB = File.stem
		Out = "Results/MotherConsistency/{SGB}.tsv".format(SGB=SGB)
		if Path(Out).exists():  continue
		Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))
		Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.motherConsistency
#SBATCH --output=logs/{SGB}.motherConsistency.o
#SBATCH --error=logs/{SGB}.motherConsistency.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=20GB
cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
		Script_n = "scripts/{SGB}_motherconsistency.sh".format(SGB = SGB)
		with open(Script_n,'w') as F: F.write(Job_setup + Command)
		job_count = check_job_count()
		while job_count >= 200:
			time.sleep(60)
			job_count = check_job_count()
		print(Script_n)
		subprocess.call("sbatch "+Script_n, shell=True)	
			
def Script_Abundance_vs_sharing():
	script ="7_RelAbundance_and_Sharing.R"
	SGBs = []
	with open("metadata/LLNEXT_metaphlan_4_CLR_transformed_fil_SGB_infants_20_07_2023.txt") as F:
		for line in F:
			l = line.rstrip().split()
			for i in l[1:]:
				SGBs.append( i.split(".")[1] )
			break
	for SGB in SGBs:
		Out = "Results/AbundanceAssociation/{S}.tsv".format(S = SGB)
		#if Path(Out).exists(): continue
		File = "Distance_tables/{S}.rds".format(S=SGB)
		if not Path(File).exists(): continue
		Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))
		Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.abundance
#SBATCH --output=logs/{SGB}.abundance.o
#SBATCH --error=logs/{SGB}.abundance.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=20GB
cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
		Script_n = "scripts/{SGB}_Abundance_vs_sharing.sh".format(SGB = SGB)
		with open(Script_n,'w') as F: F.write(Job_setup + Command)
		subprocess.call("sbatch "+Script_n, shell=True)
		print(Script_n)
			
def Script_AbundanceMother_vs_sharing():
        script ="8_RelAbundanceMother_sharing.R"
        SGBs = []
        with open("metadata/NEXT_metaphlan_4_CLR_transformed_fil_SGB_mothers_03_08_2023.txt") as F:
                for line in F:
                        l = line.rstrip().split()
                        for i in l[1:]:
                                SGBs.append( i.split(".")[1] )
                        break
        for SGB in SGBs:
                Out = "Results/AbundanceAssociationMother/{S}.tsv".format(S = SGB)
                #if Path(Out).exists(): continue
                File = "Distance_tables/{S}.rds".format(S=SGB)
                if not Path(File).exists(): continue
                Command = "ml RPlus\nRscript {script} {F}".format(script = script, F=str(File))
                Job_setup = """#!/bin/bash
#SBATCH --job-name={SGB}.abundanceMother
#SBATCH --output=logs/{SGB}.abundanceMother.o
#SBATCH --error=logs/{SGB}.abundanceMother.e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mem=20GB
cd /groups/umcg-llnext/tmp01/umcg-sandreusanchez/Strain_transmission
""".format(SGB = SGB)
                Script_n = "scripts/{SGB}_AbundanceMother_vs_sharing.sh".format(SGB = SGB)
                with open(Script_n,'w') as F: F.write(Job_setup + Command)
                subprocess.call("sbatch "+Script_n, shell=True)
                print(Script_n)



#Script1_transmission()
#Script3_mantelTest()
Script2_LinearModel_transmission()
#Secript_motherConsistency()
#Script_Abundance_vs_sharing()
#Script_AbundanceMother_vs_sharing()


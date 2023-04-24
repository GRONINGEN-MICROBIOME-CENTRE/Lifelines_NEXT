####################################################################################
# It is a script for shuffling the order of viral contigs inside of the fasta file
# so that it will take equal time to run the tool for the analyses
####################################################################################

##############################
# Uploading libraries  
##############################
import pandas as pd
import random
from pyfaidx import Fasta
import textwrap

##############################
# Uploading data  
##############################
faa = Fasta("C:\\Users\\Natal\\PycharmProjects\\iphop\\final_data\\viral_noneg405_99_der95_decontaminated.fasta")

metadata = pd.read_table("C:\\Users\\Natal\\PycharmProjects\\iphop\\metadata\\VLP_viral_contigs_metadata.txt",
                         header=0)

#########################################################################
# Selecting viral contigs for the host prediction based on metadata file 
#########################################################################
included = metadata["V1"].tolist()
random.shuffle(included)

###################################################################
# Shuffling viral contigs and saving the result in the .fasta file
###################################################################
with open("viral_noneg405_99_der95_decontaminated_filtered_shuffled.fasta", "w") as outfile:
    for i in included:
        outfile.write(f'>{i}\n')
        outfile.write("\n".join(textwrap.wrap(str(faa[i]), 60)))
        outfile.write("\n")
print("success")

######################################
# Additional check step for the result
######################################

# File check (moved obtained out files to 'results' folder)
# faa = Fasta("C:\\Users\\Natal\\PycharmProjects\\iphop\\Scripts\\viral_noneg405_99_der95_decontaminated_filtered_shuffled.fasta")
# for i in faa.keys():
#     if int(str(i).split("_")[7]) != len(faa[i]):
#         print(i)
#         print(len(faa[i]))
# print(len(faa.keys()))

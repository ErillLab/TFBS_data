import csv


with open("TFSet.txt") as f:
    lines = [line for line in csv.reader(f,delimiter='\t') if not line[0].startswith("#")]
gene_names_from_tf_dict = {line[1]:line[2].split(", ") for line in lines}

def gene_name_from_tf(tf):
    return gene_names_from_tf_dict[tf]
        



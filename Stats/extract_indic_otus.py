import pandas as pd

max_otus = 5
with open("../Data/otus_R1.fasta", "r") as ifasta:
    otu_dict = {}
    line = ifasta.readline()
    while line != "":
        name = line.strip().strip(">")
        seq = ""
        line = ifasta.readline()
        while line != ""  and ">" not in line:
            seq += line
            line = ifasta.readline()
        otu_dict[name] = seq
constax_dict = {}
with open("../Data/constax_taxonomy.txt", "r") as ifile:
    line = ifile.readline()
    while line != "":
        spl = line.strip().split("\t")
        print(spl)
        if len(spl) == 8:
            tax = spl[7].replace(" ", "_")
        else:
            tax = F"{spl[-1]}_sp"
        constax_dict[spl[0]] = tax
        line = ifile.readline()

indic_dat = pd.read_csv("OTU_indicspecies_fdr.csv")
indic_sig = indic_dat[indic_dat["p.value"] <= 0.05]
indic_sig.sort_values(["s.Endo", "s.Epi", "s.Lit", "s.Soil", "p.value"], inplace = True)
categ_dict = {}
fa_buf = ""
csv_buf = ",OTU," + ",".join(indic_dat.columns) + ",taxonomic_id\n"
for i in range(len(indic_sig)):
    id = tuple(indic_sig.iloc[i, 0:4])
    if id not in categ_dict.keys():
        categ_dict[id] = 1
    elif categ_dict[id] < max_otus:
        categ_dict[id] += 1
    else:
        continue
    name = indic_sig.index[i]
    fa_buf += F">{name}\n{otu_dict[name]}"
    csv_buf += F"{i},{name},{str(list(indic_sig.iloc[i]))[1:-1].replace(', ', ',')},{constax_dict[name]}\n"
with open("indic_species_5_per.csv", "w") as ocsv:
    ocsv.write(csv_buf)
with open("indic_species_5_per.fasta", "w") as ofasta:
    ofasta.write(fa_buf)

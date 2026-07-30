# Merging every samples files with ARG + species in 1 single file

import csv
import glob
from pathlib import Path

output_file = "all_samples_merged.csv"
first_file = True

with open(output_file, "w", newline="") as out:
    writer = None

    # Parcourt tous les fichiers merged.csv
    for merged_file in glob.glob("*_merged.csv"):

        sample = Path(merged_file).name.replace("_merged.csv", "")
        print(f"Ajout du sample : {sample}")

        with open(merged_file, newline="") as f:
            reader = csv.DictReader(f)

            # Initialisation de l'en-tête une seule fois
            if first_file:
                fieldnames = ["sample"] + reader.fieldnames
                writer = csv.DictWriter(out, fieldnames=fieldnames)
                writer.writeheader()
                first_file = False

            # Écriture des lignes avec ajout du sample
            for row in reader:
                row_with_sample = {"sample": sample}
                row_with_sample.update(row)
                writer.writerow(row_with_sample)

print(f"Final joined file : {output_file}")

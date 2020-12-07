#Lectures des fichiers .csv, extraction du nombre total des échantillons (hors contrôles)

total_sample = 0
total_retest = 0

for (file in list.files(pattern = ".csv")) {
  data <- read.csv(file = file)
  for (id in data[,1]) {
    if (substr(id,1,4) == "106C") {
      total_sample = total_sample + 1
    }
  }
  for (retest in data[,"IsRetest"]) {
    if (retest == "True") {
      total_retest = total_retest + 1
    }
  }
}

print(total_sample)
print(total_retest)

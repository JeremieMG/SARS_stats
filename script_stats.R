read.csv("Essai stat rerun.csv") -> data

#Partie 1: Calcul du nombre de changements de résultats

inc_to_inc = 0
inc_to_pos = 0
inc_to_neg = 0
pos_to_inc = 0
pos_to_pos = 0
pos_to_neg = 0
neg_to_inc = 0
neg_to_pos = 0
neg_to_neg = 0
val_to_inc = 0
val_to_pos = 0
val_to_neg = 0


for (id in data[,"Sample.ID"]) {
  if (data[which(data["Sample.ID"]==id),"Run.1"]=="Inconclusive") {
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Inconclusive") {
      inc_to_inc = inc_to_inc + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Detected") {
      inc_to_pos = inc_to_pos + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Not Detected") {
      inc_to_neg = inc_to_neg + 1
    }
  }
  if (data[which(data["Sample.ID"]==id),"Run.1"]=="Detected") {
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Inconclusive") {
      pos_to_inc = pos_to_inc + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Detected") {
      pos_to_pos = pos_to_pos + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Not Detected") {
      pos_to_neg = pos_to_neg + 1
    }
  }
  if (data[which(data["Sample.ID"]==id),"Run.1"]=="Not Detected") {
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Inconclusive") {
      neg_to_inc = neg_to_inc + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Detected") {
      neg_to_pos = neg_to_pos + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Not Detected") {
      neg_to_neg = neg_to_neg + 1
    }
  }
  else {
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Inconclusive") {
      val_to_inc = val_to_inc + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Detected") {
      val_to_pos = val_to_pos + 1
    }
    if (data[which(data["Sample.ID"]==id),"Run.2"]=="Not Detected") {
      val_to_neg = val_to_neg + 1
    }
  }
}

m1 <- matrix(nrow = 4, ncol = 3, 
       data = c(inc_to_inc, inc_to_pos, inc_to_neg, pos_to_inc, pos_to_pos, pos_to_neg, neg_to_inc, neg_to_pos, neg_to_neg, val_to_inc, val_to_pos, val_to_neg),
       dimnames = list(c("Inconclusif","Positif","Négatif","Invalide"),c("devient: Inconclusif","Positif","Négatif")),
       byrow = "TRUE")

#Partie 2: Calculs en fonction du Ct de N

which(is.na(data[,"N"]) != "TRUE" & is.na(data[,"MS2"]) != "TRUE" & is.na(data[,"ORF1ab"]) == "TRUE" & is.na(data[,"S"]) == "TRUE") -> subset

ns_inc = 0
ns_pos = 0
ns_neg = 0
ni_inc = 0
ni_pos = 0
ni_neg = 0

for (position in subset) {
  if (data[position,"N"] >= 32) {
    if (data[position,"Run.2"] == "Inconclusive") {
      ns_inc = ns_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      ns_pos = ns_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      ns_neg = ns_neg + 1 
    }
  }
  else {
    if (data[position,"Run.2"] == "Inconclusive") {
      ni_inc = ni_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      ni_pos = ni_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      ni_neg = ni_neg + 1 
    }
  }
}

m2 <-matrix(nrow = 2, ncol = 3, data = c(ns_inc/(ns_inc+ns_pos+ns_neg)*100,ns_pos/(ns_inc+ns_pos+ns_neg)*100,ns_neg/(ns_inc+ns_pos+ns_neg)*100,ni_inc/(ni_inc+ni_pos+ni_neg)*100,ni_pos/(ni_inc+ni_pos+ni_neg)*100,ni_neg/(ni_inc+ni_pos+ni_neg)*100),
       dimnames = list(c("N >=32","N <32"),c("Inc (%)","Pos (%)","Neg (%)")),
       byrow = "TRUE")

#Calculs en fonction de ORF1ab

which(is.na(data[,"N"]) == "TRUE" & is.na(data[,"MS2"]) != "TRUE" & is.na(data[,"ORF1ab"]) != "TRUE" & is.na(data[,"S"]) == "TRUE") -> subset

os_inc = 0
os_pos = 0
os_neg = 0
oi_inc = 0
oi_pos = 0
oi_neg = 0

for (position in subset) {
  if (data[position,"ORF1ab"] >= 32) {
    if (data[position,"Run.2"] == "Inconclusive") {
      os_inc = os_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      os_pos = os_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      os_neg = os_neg + 1 
    }
  }
  else {
    if (data[position,"Run.2"] == "Inconclusive") {
      oi_inc = oi_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      oi_pos = oi_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      oi_neg = oi_neg + 1 
    }
  }
}

m3 <-matrix(nrow = 2, ncol = 3, data = c(os_inc/(os_inc+os_pos+os_neg)*100,os_pos/(os_inc+os_pos+os_neg)*100,os_neg/(os_inc+os_pos+os_neg)*100,oi_inc/(oi_inc+oi_pos+oi_neg)*100,oi_pos/(oi_inc+oi_pos+oi_neg)*100,oi_neg/(oi_inc+oi_pos+oi_neg)*100),
            dimnames = list(c("ORF >=32","ORF <32"),c("Inc","Pos","Neg")),
            byrow = "TRUE")

#Calculs en fonction de S

which(is.na(data[,"N"]) == "TRUE" & is.na(data[,"MS2"]) != "TRUE" & is.na(data[,"ORF1ab"]) == "TRUE" & is.na(data[,"S"]) != "TRUE") -> subset

ss_inc = 0
ss_pos = 0
ss_neg = 0
si_inc = 0
si_pos = 0
si_neg = 0

for (position in subset) {
  if (data[position,"S"] >= 32) {
    if (data[position,"Run.2"] == "Inconclusive") {
      ss_inc = ss_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      ss_pos = ss_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      ss_neg = ss_neg + 1 
    }
  }
  else {
    if (data[position,"Run.2"] == "Inconclusive") {
      si_inc = si_inc + 1 
    }
    if (data[position,"Run.2"] == "Detected") {
      si_pos = si_pos + 1 
    }
    if (data[position,"Run.2"] == "Not detected") {
      si_neg = si_neg + 1 
    }
  }
}

m4 <-matrix(nrow = 2, ncol = 3, data = c(ss_inc/(ss_inc+ss_pos+ss_neg)*100,ss_pos/(ss_inc+ss_pos+ss_neg)*100,ss_neg/(ss_inc+ss_pos+ss_neg)*100,si_inc/(si_inc+si_pos+si_neg)*100,si_pos/(si_inc+si_pos+si_neg)*100,si_neg/(si_inc+si_pos+si_neg)*100),
            dimnames = list(c("S >=32","S <32"),c("Inc","Pos","Neg")),
            byrow = "TRUE")

print(m1)
print(m2)
print(m3)
print(m4)
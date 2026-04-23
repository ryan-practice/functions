library(usethis)
library(data.table)
library(functions)

traj <- fread(wsl_ify("C:/Users/rms2jg/Documents/early_onboard/TRAJ.csv"), sep = ",")
traj$Search <- ""
traj$Output <- ""
traj$n_gxg <- sapply(traj$`AA Seq`, function(s) {
                              m <- gregexpr("G.G", s)[[1]]
                              sum(m > 0)
                            })
traj$aa_of_traj_1 <- ""
traj$aa_of_traj_2 <- ""
traj[which(traj$TRAJ == "TRAJ16"), "aa_of_traj_1"] <- "XFSDGQKLL"
traj[which(traj$TRAJ == "TRAJ16"), "aa_of_traj_2"] <- "FARGTMLKVDL"
traj[which(traj$TRAJ == "TRAJ28"), "aa_of_traj_1"] <- "XYSGAGSYQLTF"
traj[which(traj$TRAJ == "TRAJ28"), "aa_of_traj_2"] <- "GKGTKLSVIP"
traj[which(traj$TRAJ == "TRAJ45"), "aa_of_traj_1"] <- "XYSGGGADGLTF"
traj[which(traj$TRAJ == "TRAJ45"), "aa_of_traj_2"] <- "GKGTHLIIQP"
trbj <- fread(wsl_ify("C:/Users/rms2jg/Documents/early_onboard/TRBJ.csv"), sep = ",")
trbj$Search <- ""
trbj$Output <- ""
trbj$n_gxg <- sapply(trbj$`AA Seq`, function(s) {
  m <- gregexpr("G.G", s)[[1]]
  sum(m > 0)
})
trav <- fread(wsl_ify("C:/Users/rms2jg/Documents/early_onboard/TRAV.csv"), sep = ",")
trav$Location <- ""
trav$`Plasmid name` <- ""
trav$`Full Plasmid name` <- ""
trbv <- fread(wsl_ify("C:/Users/rms2jg/Documents/early_onboard/TRBV.csv"), sep = ",")
trbv$Location <- ""
trbv$`Plasmid name` <- ""
trbv$`Full Plasmid name` <- ""

usethis::use_data(traj, internal = TRUE, overwrite = TRUE)
usethis::use_data(trbj, internal = TRUE, overwrite = TRUE)
usethis::use_data(trav, internal = TRUE, overwrite = TRUE)
usethis::use_data(trbv, internal = TRUE, overwrite = TRUE)

template_fasta <- fread(wsl_ify("C:/Users/rms2jg/Documents/early_onboard/eblock-template-cdr3s.fasta"), sep = "")
template_fasta <- paste(template_fasta$`>eBlock template-CDR3s`, collapse = "")

usethis::use_data(template_fasta, overwrite = TRUE)

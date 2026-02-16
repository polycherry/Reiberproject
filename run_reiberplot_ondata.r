wd = "/Volumes/lab-schmackk/home/shared/rawData/009_antibodyAssays/AbAs001_IRIS/AbAs001_06_CSF_reiberanalysis/Reiber project"
source(file.path(wd, "igg_reiberplot.r"))


# Example patient data
dfa <- data.frame(
  patient_id = c("101a", "103a", "104a", "105a", "107a", "108a", "109a", "112a", "113a", "114a", "117a", "120a", "121a"),
  csf_igg = c(25.3, 22.2, 50.9, 24.9, 56.9, 16.7, 23.3, 19.9, 36.9, 17, 22.3, 35, 17.7),
  s_igg = c(13, 9.16, 17.5, 10, 13.1, 10.4, 16.6, 9.83, 18, 13.5, 12.4, 15.9, 17.3),
  csf_alb = c(179, 213, 260, 196, 373, 149, 130, 176, 203, 117, 170, 204, 91.6),
  s_alb = c(45.3, 42.7, 41.9, 46.8, 46.1, 43.3, 44.4, 43.8, 49.5, 41, 43.6, 46.2, 43.8),
  age = c(33, 31, 27, 23, 28, 26, 24, 21, 18, 21, 46, 22, 18)
)

dfa$csf_igg = dfa$csf_igg/1000
dfa$csf_alb = dfa$csf_alb/1000
# Plot
plot_reibergram(dfa)

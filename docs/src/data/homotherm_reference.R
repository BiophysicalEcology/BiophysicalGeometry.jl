# Reference geometry of the default NicheMapR human, for the tutorial "A human: comparison with NicheMapR".
#
# Run once with NicheMapR installed; the docs read the CSV files this writes.
#   Rscript homotherm_reference.R

library(NicheMapR)

MASS <- 70
parts <- c("head", "trunk", "arm", "leg")
NPARTs <- c(1, 1, 2, 2)
SHAPEs <- c(4, 1, 1, 1) # ellipsoid head, cylinder trunk, arms and legs
MASSFRACs <- c(0.07609801, 0.50069348, 0.04932963, 0.16227462)
DENSITYs <- rep(1050, 4)
SHAPE_Bs <- c(1.6, 1.9, 12, 7.0)
FATPCTs <- c(5, 36, 10, 23) * 0.7
PJOINs <- c(0.02666667, 0.08088128, 0.02000000, 0.03333333)
INSDEPDs <- c(1e-02, rep(6e-03, 3))
INSDEPVs <- c(1e-09, rep(6e-03, 3))
INSDEPs <- (INSDEPDs + INSDEPVs) / 2 # the mean depth HomoTherm uses for geometry
MASSs <- MASS * MASSFRACs

GEOM.lab <- c("VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH", "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP",
              "GMASS", "AREASKIN", "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK", "CONVAR", "R1", "R2")
geom <- function(i, ZFUR, ORIENT = 0, ZEN = 0) {
  g <- GEOM_ENDO(MASSs[i], DENSITYs[i], DENSITYs[i], FATPCTs[i], SHAPEs[i], ZFUR, 1, SHAPE_Bs[i], SHAPE_Bs[i],
                 0, 0, PJOINs[i], 0, ORIENT, ZEN)
  setNames(as.numeric(g), GEOM.lab)
}

# per-part geometry, direct from GEOM_ENDO
g <- t(sapply(1:4, function(i) geom(i, INSDEPs[i])))
skin_length <- ifelse(SHAPEs == 4, 2 * g[, "ASEMAJ"], 2 * SHAPE_Bs * g[, "R1"])

# per-part geometry as reported by the full model
out <- HomoTherm()
morph <- rbind(out$head.morph, out$trunk.morph, out$arm.morph, out$leg.morph)

part_table <- data.frame(
  part = parts, n = NPARTs, mass_kg = MASSs, density_kg_m3 = DENSITYs, shape_b = SHAPE_Bs, fat_pct = FATPCTs,
  insulation_depth_m = INSDEPs, pjoin = PJOINs,
  volume_m3 = g[, "VOL"], flesh_volume_m3 = g[, "FLSHVL"], fat_thickness_m = g[, "FATTHK"],
  skin_radius_m = g[, "R1"], insulation_radius_m = g[, "R2"], semi_major_m = g[, "ASEMAJ"],
  skin_length_m = skin_length, skin_area_m2 = g[, "AREASKIN"], total_area_m2 = g[, "ATOT"],
  homotherm_area_m2 = morph[, "AREA"], homotherm_skin_area_m2 = morph[, "AREA_SKIN"],
  homotherm_join_area_m2 = morph[, "AREA_JOIN"])
write.csv(part_table, "homotherm_parts.csv", row.names = FALSE)

# whole-body values
height <- 2 * g[1, "ASEMAJ"] + skin_length[2] + skin_length[4] + INSDEPDs[4] # as returned by plot_human
totals <- data.frame(
  quantity = c("area_net_m2", "area_gross_m2", "skin_area_gross_m2", "join_area_m2", "height_m", "dubois_area_m2"),
  value = c(out$balance["AREA"], sum(NPARTs * g[, "ATOT"]), sum(NPARTs * g[, "AREASKIN"]),
            sum(NPARTs * morph[, "AREA_JOIN"]), height, 0.00718 * MASS^0.425 * (100 * height)^0.725))
write.csv(totals, "homotherm_totals.csv", row.names = FALSE)

# silhouette area against zenith angle, as HomoTherm calls it
ZENs <- seq(0, 90, 5)
sil <- NicheMapR:::human_silhouette_area(ZENs = ZENs, MASSs = MASSs, DENSITYs = DENSITYs, SHAPE_Bs = SHAPE_Bs,
                                         FATPCTs = FATPCTs, PJOINs = PJOINs, plot.sil = FALSE)
write.csv(data.frame(zenith_deg = ZENs, parts_sum_m2 = sil$sil.HomoTherm, underwood_ward_m2 = sil$sil.underwood),
          "homotherm_silhouette.csv", row.names = FALSE)

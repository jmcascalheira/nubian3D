
## EXTENDED METHODS STEP 8: R SCRIPT FOR GEOMETRIC MORPHOMETRIC ANALYSES

# Ensure all require libraries are loaded
library(geomorph)
library(Rvcg)
library(ggplot2)
library(car)

# If not already imported in Step 7:
## -- CORE DATA --
# Core surfaces
load("analysis/data/derived_data/allpatched.RData")

# Core outlines
load("analysis/data/derived_data/alloutline.RData")

# Pref scar outlines
load("analysis/data/derived_data/allprefs.RData")

## -- PRODUCT DATA -- (Not used in Step 7)
# NK end products
load("analysis/data/derived_data/prod_patched.Rdata")
load("analysis/data/derived_data/prod_ventr.Rdata")

# Import attribute files
cores <- read.csv("analysis/data/derived_data/File_list_cores.csv")
cores <- dplyr::arrange(cores, File_name)

prods <- read.csv("analysis/data/derived_data/File_list_prods.csv")
prods <- dplyr::arrange(prods, File_name)

###---------------------------------------------------------------------------

## STEP 8.1 Conduct Generalised Procrustes Analysis to scale, rotate and align shapes (gpagen in geomorph)
# GPA for cores already run in Step 7:
patched_gpa <- gpagen(allpatched,PrinAxes=TRUE)
outline_gpa <- gpagen(alloutline,PrinAxes=TRUE)
pref_gpa <- gpagen(allprefs,PrinAxes=TRUE)

# Run GPA for products
prod_gpa <- gpagen(prod_patched, PrinAxes=TRUE)
prodv_gpa <- gpagen(prod_ventr, PrinAxes=TRUE)

###---------------------------------------------------------------------------

## STEP 8.2 Identify and exclude any outlier specimens based on Procrustes distance from the mean shape
## -- CORE SURFACES --
# Calculate Procrustes distances to the mean shape
meansh <- mshape(patched_gpa$coords)
proc_dist <- apply(patched_gpa$coords, 3, function(x) sqrt(sum((x - meansh)^2)))

# Calculate the mean and standard deviation of the Procrustes distances
mean_proc <- mean(proc_dist)
sd_proc <- sd(proc_dist)

# Define the threshold for outliers within 3 sd of the mean (i.e. 99.7% of the data)
threshold_upper <- mean_proc + 3 * sd_proc
threshold_lower <- mean_proc - 3 * sd_proc

# Identify outliers
outliers <- proc_dist < threshold_lower | proc_dist > threshold_upper

# Display outlier specimens
outlier_specs <- proc_dist[outliers]
print(outlier_specs)

# Exclude outliers from the patched core dataset and rerun GPA on new dataset
outlier_specs <- c(70, 91, 72, 119) # identify corresponding specimen numbers from corelist
all_specs <- 1:166
keep <- setdiff(all_specs, outlier_specs)
allpatched <- allpatched[, , keep]
patched_gpa <- gpagen(allpatched,PrinAxes=TRUE) # check the new patched_gpa has 162 specimens


## -- PREF SCAR OUTLINES --
# Calculate Procrustes distances to the mean shape
meansh <- mshape(pref_gpa$coords)
proc_dist <- apply(pref_gpa$coords, 3, function(x) sqrt(sum((x - meansh)^2)))

# Calculate the mean and standard deviation of the Procrustes distances
mean_proc <- mean(proc_dist)
sd_proc <- sd(proc_dist)

# Define the threshold for outliers within 3 sd of the mean
threshold_upper <- mean_proc + 3 * sd_proc
threshold_lower <- mean_proc - 3 * sd_proc

# Identify outliers
outliers <- proc_dist < threshold_lower | proc_dist > threshold_upper

# Display outlier specimens
outlier_specs <- proc_dist[outliers]
print(outlier_specs)

# Exclude outliers from the pref scar dataset and rerun GPA on new dataset
outlier_specs <- c(45) # identify corresponding specimen numbers from corelist
all_specs <- 1:166
keep <- setdiff(all_specs, outlier_specs)
allprefs <- allprefs[, , keep]
pref_gpa <- gpagen(allprefs,PrinAxes=TRUE) # check the new pref_gpa has 165 specimens

# Exclude corresponding outliers from the core outline dataset (needed for 2BPLS between shapes on the same specimens) and rerun GPA on new dataset
alloutline <- alloutline[, , keep]
outline_gpa <- gpagen(alloutline,PrinAxes=TRUE) # check the new pref_gpa has 165 specimens


## -- PRODUCT SURFACES --
# Calculate Procrustes distances to the mean shape
meansh <- mshape(prod_gpa$coords)
proc_dist <- apply(prod_gpa$coords, 3, function(x) sqrt(sum((x - meansh)^2)))

# Calculate the mean and standard deviation of the Procrustes distances
mean_proc <- mean(proc_dist)
sd_proc <- sd(proc_dist)

# Define the threshold for outliers within 3 sd of the mean
threshold_upper <- mean_proc + 3 * sd_proc
threshold_lower <- mean_proc - 3 * sd_proc

# Identify outliers
outliers <- proc_dist < threshold_lower | proc_dist > threshold_upper

# Display outlier specimens
outlier_specs <- proc_dist[outliers]
print(outlier_specs)

# Exclude outliers from the patched product dataset and rerun GPA on new dataset
outlier_specs <- c(8, 28, 96) # identify corresponding specimen numbers from prodlist
all_specs <- 1:178
keep <- setdiff(all_specs, outlier_specs)
prod_patched <- prod_patched[, , keep]
prod_gpa <- gpagen(prod_patched,PrinAxes=TRUE) # check the new prod_gpa has 175 specimens

# Exclude the corresponding outliers from the product ventral dataset and rerun GPA on new dataset
prod_ventr <- prod_ventr[, , keep]
prodv_gpa <- gpagen(prod_ventr,PrinAxes=TRUE) # check the new prodv_gpa has 175 specimens

###---------------------------------------------------------------------------

## STEP 8.3 Conduct PCAs (gm.prcomp in geomorph)
## -- CORE SURFACES --
# Conduct PCA and extract PC scores to new dataframe
PCA_allpatch <- gm.prcomp(patched_gpa$coords)
PCscores_patch <- PCA_allpatch$x # x represents the PC scores
coord_allpatch <- as.data.frame(PCscores_patch)

# Specify outlier specimen names and exclude these specimens from the attribute dataset
outlier_specs <- c("CoreME78.737b_LR", "CoreME78.959_LR", "CoreME78.747.78_LR", "TH.571-110_LR")
cores_filtered <- cores[!cores$File_name %in% outlier_specs, ]

# Append filtered core attribute data to the coord data
coord_allpatch <- cbind(coord_allpatch, cores_filtered)

# Generate variance tables that show the percentage of variation captured by each component
eigenvalues <- PCA_allpatch$d # d represents the eigenvalues
pc_numbers <- 1:length(eigenvalues)
perc_variance <- (eigenvalues / sum(eigenvalues)) * 100
cum_perc <- cumsum(perc_variance)
perc_table_patch <- data.frame(PC = 1:length(eigenvalues), Eigenvalue = eigenvalues, PercVariance = perc_variance, CumulativePerc = cum_perc)


## -- PREF SCAR OUTLINE --
# Conduct PCA and extract PC scores to new dataframe
PCA_allprefs <- gm.prcomp(pref_gpa$coords)
PCscores_prefs <- PCA_allprefs$x
coord_allprefs <- as.data.frame(PCscores_prefs)

# Specify outliers and exclude these specimens from the attribute dataset
outlier_specs <- c("CoreME78.627Ap_LR")
cores_filtered <- cores[!cores$File_name %in% outlier_specs, ]

# Append filtered core attribute data to the coord data
coord_allprefs <- cbind(coord_allprefs, cores_filtered)

# Generate variance tables that show the percentage of variation captured by each component
eigenvalues <- PCA_allprefs$d
pc_numbers <- 1:length(eigenvalues)
perc_variance <- (eigenvalues / sum(eigenvalues)) * 100
cum_perc <- cumsum(perc_variance)
perc_table_pref <- data.frame(PC = 1:length(eigenvalues), Eigenvalue = eigenvalues, PercVariance = perc_variance, CumPercentage = cum_perc)


## -- PRODUCT SURFACES --
# Conduct PCA and extract PC scores to new dataframe
PCA_prods <- gm.prcomp(prod_gpa$coords)
PCscores_prods <- PCA_prods$x
coord_prods <- as.data.frame(PCscores_prods)

# Specify outliers and exclude these specimens from the attribute dataset
outlier_specs <- c("PointME78.1049c_LR", "PointME78.369.8_LR", "PointME78.718p_LR")
prods_filtered <- prods[!prods$File_name %in% outlier_specs, ]

# Append filtered product attribute data to the coord data
coord_prods <- cbind(coord_prods, prods_filtered)

# Generate variance tables that show the percentage of variation captured by each component
eigenvalues <- PCA_prods$d
pc_numbers <- 1:length(eigenvalues)
perc_variance <- (eigenvalues / sum(eigenvalues)) * 100
cum_perc <- cumsum(perc_variance)
perc_table_prods <- data.frame(PC = 1:length(eigenvalues), Eigenvalue = eigenvalues, PercVariance = perc_variance, CumPercentage = cum_perc)

###---------------------------------------------------------------------------

## STEP 8.4 Plot PCAs. Repeat for products and showing other variables (e.g. core preparation) as aesthetics
# If not done in Step 7, set colour/shape aesthetics for the NK and TH samples for the plots
colsArea <- c("NK" = "royalblue", "TH" = "indianred1")
shapesArea <- c("NK" = 16, "TH" = 17)

## -- CORE SURFACES --
## Generate PCA plots (here with aesthetics set for area/region)
#PC1 and 2 (Figure 5)
coord_allpatch %>%
  ggplot(aes(x = Comp1, y = Comp2, shape = Area)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Area)) +
  scale_color_manual(values = colsArea) +
  scale_fill_manual(values = colsArea) +
  stat_ellipse(geom = "polygon", aes(fill = Area), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(shape = "Region", color = "Region", fill = "Region") +
  labs(x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("Core surfaces: PCA (PC1-PC2)")

#PC3 and 4 (Figure 6)
coord_allpatch %>%
  ggplot(aes(x = Comp3, y = Comp4, shape = Area)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Area)) +
  scale_color_manual(values = colsArea) +
  scale_fill_manual(values = colsArea) +
  stat_ellipse(geom = "polygon", aes(fill = Area), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(shape = "Region", color = "Region", fill = "Region") +
  labs(x = "PC3", y = "PC4") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("Core surfaces: PCA (PC3-PC4)")

#PC5 and 6
coord_allpatch %>%
  ggplot(aes(x = Comp5, y = Comp6, shape = Area)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Area)) +
  scale_color_manual(values = colsArea) +
  scale_fill_manual(values = colsArea) +
  stat_ellipse(geom = "polygon", aes(fill = Area), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(shape = "Region", color = "Region", fill = "Region") +
  labs(x = "PC5", y = "PC6") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("Core surfaces: PCA (PC5-PC6)")

# Generate PCA plots with aesthetics for Type and Assemblage (see supplementary script 8A)

###---------------------------------------------------------------------------

## STEP 8.5 Additional visualisation steps to represent shape variation of core surfaces
# Find mean, minimum and maximum shapes (mshape and PlotReftoTarget visualisation in geomorph)

## -- CORE SURFACES --
# Find the theoretical mean shape
msh_patch <- mshape(patched_gpa$coords)

# Create an outline by linking consecutive slms of the shape
#define.links(msh_patch, ptsize = 3, links = NULL) # can define links between slms manually or input as a .csv
core_links <- read.csv("analysis/data/derived_data/Core_links.csv")

# Identify the specimen that most closely matches the theoretical mean shape
findMeanSpec(patched_gpa$coords)
#CoreME78.789.41_LR
#Specimen 83


# Define a function to generate and save min./max. shape snapshots for a given PC
GeneratePCShapes <- function(pc_num) {
  # Mean shape (grey) and minimum shape (blue)
  plotRefToTarget(msh_patch, PCA_allpatch$shapes[[paste0("shapes.comp", pc_num)]]$min,
                  method = "points", links = core_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "blue", tar.pt.size = 0.7,
                  tar.link.col = "blue", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_min_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_min.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_min_dist.png"), fmt = 'png')

  # Mean shape (grey) and maximum shape (red)
  plotRefToTarget(msh_patch, PCA_allpatch$shapes[[paste0("shapes.comp", pc_num)]]$max,
                  method = "points", links = core_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "red", tar.pt.size = 0.7,
                  tar.link.col = "red", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_max_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_max.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("PC", pc_num, "_max_dist.png"), fmt = 'png')
}

# Loop through PCs 1 to 6 (or desired number of PCs)
for (pc in 1:6) {
  GeneratePCShapes(pc)
}


# Identify the specimens corresponding to the minimum and maximum shapes
# Create a new dataframe to hold the extracted specimen information
Specs <- data.frame(PC = character(), MinSp = character(), MaxSp = character(), stringsAsFactors = FALSE)

# Identify specimen names in coord_allpatch dataframe
coord_allpatch <- data.frame(
  Specimen = rownames(coord_allpatch),
  coord_allpatch,
  stringsAsFactors = FALSE)

# Extract specimen names that correspond to min. and max. value for first 6 PCs
for (i in 1:6) {
  pc_name <- paste0("Comp", i)
  min_specimen <- coord_allpatch$Specimen[which.min(coord_allpatch[[pc_name]])]
  max_specimen <- coord_allpatch$Specimen[which.max(coord_allpatch[[pc_name]])]
  Specs <- rbind(Specs, data.frame(PC = pc_name, MinSp = min_specimen, MaxSp = max_specimen, stringsAsFactors = FALSE))
}
print(Specs)


# Refer to specimens in the table to identify the correct mesh and outline slms (get specimen number from corelist) and output views for min. and max. specimens
# Example for PC1 minimum specimen
corePC1min <- vcgImport(file.path("analysis", "data", "raw_data", "core_data", "TH.571-130_LR.ply"))
shade3d(corePC1min, color = "gray")
spheres3d(alloutline[1:2,,130], color = "black", radius = 0.7)
spheres3d(alloutline[3,,130], color = "grey40", radius = 0.7)
spheres3d(alloutline[4:39,,130], color = "blue", radius = 0.5)
#spheres3d(allpatched[,,130], color = "blue", radius = 0.4) # view surface slm configuration

# Manually reorient and save outputs
#rgl.snapshot('PC1sp_min.png', fmt = 'png')
#rgl.snapshot('PC1sp_min_profile.png', fmt = 'png')
#rgl.snapshot('PC1sp_min_dist.png', fmt = 'png')


# Example for PC1 maximum specimen
corePC1max <- vcgImport(file.path("analysis", "data", "raw_data", "core_data", "CoreME78.987c_LR.ply"))
shade3d(corePC1max, color = "gray")
spheres3d(alloutline[1:2,,93], color = "black", radius = 0.7)
spheres3d(alloutline[3,,93], color = "grey40", radius = 0.7)
spheres3d(alloutline[4:39,,93], color = "red", radius = 0.5)
#spheres3d(allpatched[,,93], color = "red", radius = 0.4) # view surface slm configuration

# Manually reorient and save outputs
#rgl.snapshot('PC1sp_max.png', fmt = 'png')
#rgl.snapshot('PC1sp_max_profile.png', fmt = 'png')
#rgl.snapshot('PC1sp_max_dist.png', fmt = 'png')

## Repeat for other PCs

###---------------------------------------------------------------------------

## STEP 8.6 Additional visualisation steps to represent shape variation of core pref scars
# Find mean, minimum and maximum shapes (mshape and PlotReftoTarget visualisation in geomorph)

## -- PREF SCAR --
# Find the theoretical mean shape
msh_pref <- mshape(pref_gpa$coords)

# Create an outline by linking consecutive slms of the shape
#define.links(msh_pref, ptsize = 3, links = NULL) # can define links between slms manually or input as a .csv
pref_links <- read.csv("analysis/data/derived_data/Pref_links.csv")

# Identify the specimen that most closely matches the theoretical mean shape
findMeanSpec(pref_gpa$coords)
#CoreME80.4.42b_LR
#Specimen 106


# Define a function to generate and save min./max. shape snapshots for a given PC
GeneratePCShapes <- function(pc_num) {
  # Mean shape (grey) and minimum shape (blue)
  plotRefToTarget(msh_pref, PCA_allprefs$shapes[[paste0("shapes.comp", pc_num)]]$min,
                  method = "points", links = pref_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "blue", tar.pt.size = 0.7,
                  tar.link.col = "blue", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_min_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_min.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_min_dist.png"), fmt = 'png')

  # Mean shape (grey) and maximum shape (red)
  plotRefToTarget(msh_pref, PCA_allprefs$shapes[[paste0("shapes.comp", pc_num)]]$max,
                  method = "points", links = pref_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "red", tar.pt.size = 0.7,
                  tar.link.col = "red", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_max_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_max.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Pref_PC", pc_num, "_max_dist.png"), fmt = 'png')
}

# Loop through PCs 1 to 6 (or desired number of PCs)
for (pc in 1:6) {
  GeneratePCShapes(pc)
}


# Identify the specimens corresponding to the minimum and maximum shapes
# Create a new dataframe to hold the extracted specimen information
PrefSpecs <- data.frame(PC = character(), MinSp = character(), MaxSp = character(), stringsAsFactors = FALSE)

# Identify specimen names in coord_allpatch dataframe
coord_allprefs <- data.frame(
  Specimen = rownames(coord_allprefs),
  coord_allprefs,
  stringsAsFactors = FALSE)

# Extract specimen names that correspond to min. and max. value for first 6 PCs
for (i in 1:6) {
  pc_name <- paste0("Comp", i)
  min_specimen <- coord_allprefs$Specimen[which.min(coord_allprefs[[pc_name]])]
  max_specimen <- coord_allprefs$Specimen[which.max(coord_allprefs[[pc_name]])]
  PrefSpecs <- rbind(PrefSpecs, data.frame(PC = pc_name, MinSp = min_specimen, MaxSp = max_specimen, stringsAsFactors = FALSE))
}
print(PrefSpecs)


# Refer to specimens in the table to identify the correct mesh and outline slms (get specimen number from corelist) and output views for min. and max. specimens
# Example for PC1 minimum specimen
prefPC1min <- vcgImport(file.path("analysis", "data", "raw_data", "core_data", "TH.584-Gu5_LR.ply"))
shade3d(prefPC1min, color = "gray")
spheres3d(allprefs[31,,163], color = "black", radius = 0.7)
spheres3d(allprefs[37,,163], color = "black", radius = 0.7)
spheres3d(allprefs[1:30,,163], color = "blue", radius = 0.5)
spheres3d(allprefs[32:36,,163], color = "blue", radius = 0.5)

# Manually reorient and save outputs
#rgl.snapshot('PrefPC1sp_min.png', fmt = 'png')
#rgl.snapshot('PrefPC1sp_min_profile.png', fmt = 'png')
#rgl.snapshot('PrefPC1sp_min_dist.png', fmt = 'png')


# Example for PC1 maximum specimen
prefPC1max <- vcgImport(file.path("analysis", "data", "raw_data", "core_data", "CoreME78.787e_LR.ply"))
shade3d(prefPC1max, color = "gray")
spheres3d(allprefs[31,,82], color = "black", radius = 0.7)
spheres3d(allprefs[37,,82], color = "black", radius = 0.7)
spheres3d(allprefs[1:30,,82], color = "red", radius = 0.5)
spheres3d(allprefs[32:36,,82], color = "red", radius = 0.5)

# Manually reorient and save outputs
#rgl.snapshot('PrefPC1sp_max.png', fmt = 'png')
#rgl.snapshot('PrefPC1sp_max_profile.png', fmt = 'png')
#rgl.snapshot('PrefPC1sp_max_dist.png', fmt = 'png')

## Repeat for other PCs

## Generate PCA plots for core pref scar (see supplementary script 8A)

###---------------------------------------------------------------------------

## STEP 8.7 Additional visualisation steps to represent shape variation of products
# Find mean, minimum and maximum shapes (mshape and PlotReftoTarget visualisation in geomorph)

## -- PRODUCT SURFACES --
# Find the theoretical mean shape
msh_prod <- mshape(prod_gpa$coords)

# Create an outline by linking consecutive slms of the shape
#define.links(msh_prod, ptsize = 3, links = NULL) # can define links between slms manually or input as a .csv
prod_links <- read.csv("analysis/data/derived_data/Prod_links.csv")

# Identify the specimen that most closely matches the theoretical mean shape

findMeanSpec(prod_gpa$coords)
#PointME78.578.3_LR
#Specimen 54


# Define a function to generate and save min./max. shape snapshots for a given PC
GeneratePCShapes <- function(pc_num) {
  # Mean shape (grey) and minimum shape (blue)
  plotRefToTarget(msh_prod, PCA_prods$shapes[[paste0("shapes.comp", pc_num)]]$min,
                  method = "points", links = prod_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "blue", tar.pt.size = 0.7,
                  tar.link.col = "blue", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_min_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_min.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_min_dist.png"), fmt = 'png')

  # Mean shape (grey) and maximum shape (red)
  plotRefToTarget(msh_prod, PCA_prods$shapes[[paste0("shapes.comp", pc_num)]]$max,
                  method = "points", links = prod_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "red", tar.pt.size = 0.7,
                  tar.link.col = "red", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_max_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_max.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prod_PC", pc_num, "_max_dist.png"), fmt = 'png')
}

# Loop through PCs 1 to 6 (or desired number of PCs)
for (pc in 1:6) {
  GeneratePCShapes(pc)
}


# Identify the specimens corresponding to the minimum and maximum shapes
# Create a new dataframe to hold the extracted specimen information
ProdSpecs <- data.frame(PC = character(), MinSp = character(), MaxSp = character(), stringsAsFactors = FALSE)

# Identify specimen names in coord_allpatch dataframe
coord_prods <- data.frame(
  Specimen = rownames(coord_prods),
  coord_prods,
  stringsAsFactors = FALSE)

# Extract specimen names that correspond to min. and max. value for first 6 PCs
for (i in 1:6) {
  pc_name <- paste0("Comp", i)
  min_specimen <- coord_prods$Specimen[which.min(coord_prods[[pc_name]])]
  max_specimen <- coord_prods$Specimen[which.max(coord_prods[[pc_name]])]
  ProdSpecs <- rbind(ProdSpecs, data.frame(PC = pc_name, MinSp = min_specimen, MaxSp = max_specimen, stringsAsFactors = FALSE))
}
print(ProdSpecs)

# Refer to specimens in the table to identify the correct mesh and outline slms (get specimen number from corelist) and output views for min. and max. specimens
# Example for PC1 minimum specimen
prodPC1min <- vcgImport(file.path("analysis", "data", "raw_data", "product_data","PointME78.589_LR.ply"))
shade3d(prodPC1min, color = "gray")
spheres3d(prod_ventr[31,,58], color = "black", radius = 0.7)
spheres3d(prod_ventr[37,,58], color = "black", radius = 0.7)
spheres3d(prod_ventr[1:30,,58], color = "blue", radius = 0.5)
spheres3d(prod_ventr[32:36,,58], color = "blue", radius = 0.5)

# Manually reorient and save outputs
#rgl.snapshot('ProdPC1sp_min.png', fmt = 'png')
#rgl.snapshot('ProdPC1sp_min_profile.png', fmt = 'png')
#rgl.snapshot('ProdPC1sp_min_dist.png', fmt = 'png')


# Example for PC1 maximum specimen
prodPC1max <- vcgImport(file.path("analysis", "data", "raw_data", "product_data","PointME80.4.31c_LR.ply"))
shade3d(prodPC1max, color = "gray")
spheres3d(prod_ventr[31,,166], color = "black", radius = 0.7)
spheres3d(prod_ventr[37,,166], color = "black", radius = 0.7)
spheres3d(prod_ventr[1:30,,166], color = "red", radius = 0.5)
spheres3d(prod_ventr[32:36,,166], color = "red", radius = 0.5)

# Manually reorient and save outputs
#rgl.snapshot('ProdPC1sp_max.png', fmt = 'png')
#rgl.snapshot('ProdPC1sp_max_profile.png', fmt = 'png')
#rgl.snapshot('ProdPC1sp_max_dist.png', fmt = 'png')

## Repeat for other PCs

## Generate PCA plots for NK products (see supplementary script 8A)


## -- PREF SCAR AND PROD VENTRAL --
# Combine pref scar and product ventral outline coordinates into a single array
NKprefs <- allprefs[ , , 1:114]
prod_scar <- array(c(prod_ventr, NKprefs), dim = c(37, 3, 289))
dimnames(prod_scar) <- list(NULL, NULL, c(dimnames(prod_ventr)[[3]], dimnames(NKprefs)[[3]]))

# Conduct GPA
prod_scar_gpa <- gpagen(prod_scar,PrinAxes=TRUE)

# Conduct PCA and extract PC scores to new dataframe
PCA_prod_scar <- gm.prcomp(prod_scar_gpa$coords)
PCscores_prod_scar <- PCA_prod_scar$x # x represents the PC scores
coord_prod_scar <- as.data.frame(PCscores_prod_scar)
coord_prod_scar$File_name <- rownames(coord_prod_scar)

# Import the combined attribute data (outliers are excluded) and set factors
pref_prod <- read.csv("analysis/data/derived_data/Pref_prod_list.csv")
pref_prod <- dplyr::arrange(pref_prod, File_name)
Assemblage <- factor(pref_prod$Assemblage)
Type <- factor(pref_prod$Artefact)
Scars <- factor(pref_prod$Scars)

# Join attribute data to PC scores
coord_prod_scar <- merge(coord_prod_scar, pref_prod, by = "File_name", all.x = TRUE)

# Generate variance tables that show the percentage of variation captured by each component
eigenvalues <- PCA_prod_scar$d
pc_numbers <- 1:length(eigenvalues)
perc_variance <- (eigenvalues / sum(eigenvalues)) * 100
cum_perc <- cumsum(perc_variance)
perc_table_prod_scar <- data.frame(PC = 1:length(eigenvalues), Eigenvalue = eigenvalues, PercVariance = perc_variance, CumPercentage = cum_perc)


# Identify mean shape
mshape_prodscar <- mshape(prod_scar_gpa$coords)

# Define a function to generate and save min./max. shape snapshots for a given PC
GeneratePCShapes <- function(pc_num) {
  # Mean shape (grey) and minimum shape (blue)
  plotRefToTarget(mshape_prodscar, PCA_prod_scar$shapes[[paste0("shapes.comp", pc_num)]]$min,
                  method = "points", links = pref_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "blue", tar.pt.size = 0.7,
                  tar.link.col = "blue", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_min_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_min.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_min_dist.png"), fmt = 'png')

  # Mean shape (grey) and maximum shape (red)
  plotRefToTarget(mshape_prodscar, PCA_prod_scar$shapes[[paste0("shapes.comp", pc_num)]]$max,
                  method = "points", links = pref_links,
                  gridPars = gridPar(pt.bg = "grey", pt.size = 0.5, link.lwd = 2,
                  tar.pt.bg = "red", tar.pt.size = 0.7,
                  tar.link.col = "red", tar.link.lwd = 2))

  # Save output in multiple orientations
  view3d(theta = 180, phi = -90, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_max_profile.png"), fmt = 'png')
  view3d(theta = 180, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_max.png"), fmt = 'png')
  view3d(theta = 90, phi = 0, zoom = 0.8, fov = 0, interactive = TRUE)
  rgl.snapshot(paste0("Prodscar", pc_num, "_max_dist.png"), fmt = 'png')
}

# Loop through PCs 1 to 4
for (pc in 1:4) {
  GeneratePCShapes(pc)
}

# Identify the specimens corresponding to the minimum and maximum shapes
# Create a new dataframe
PrefSpecs <- data.frame(PC = character(), MinSp = character(), MaxSp = character(), stringsAsFactors = FALSE)

# Extract specimen names that correspond to min. and max. value for first 6 PCs
for (i in 1:4) {
  pc_name <- paste0("Comp", i)
  min_specimen <- coord_prod_scar$File_name[which.min(coord_prod_scar[[pc_name]])]
  max_specimen <- coord_prod_scar$File_name[which.max(coord_prod_scar[[pc_name]])]
  PrefSpecs <- rbind(PrefSpecs, data.frame(PC = pc_name, MinSp = min_specimen, MaxSp = max_specimen, stringsAsFactors = FALSE))
}
print(PrefSpecs)

## Generate PCA plots for NK pref scars and products (see supplementary script 8A)

###---------------------------------------------------------------------------

## STEP 8.8 Perform ANOVA (Type III) and Tukey post-hoc tests to examine relationships within groups of variables
## -- CORE SURFACES --
# Define list of PCs and categorical variables
pcs <- paste0("Comp", 1:6)
variables <- c("Area", "Type", "Assemblage", "Preparation")

## Create function to perform ANOVA and Tukey post-hoc tests on each PC (1-6)
AnovaTukey <- function(var) {
  for (pc in pcs) {
    formula <- as.formula(paste(pc, "~", var))
    aov_model <- aov(formula, data = coord_allpatch)
    anova_result <- Anova(aov_model, type = "III")
    print(anova_result)
    if (length(unique(coord_allpatch[[var]])) > 2) {
      tukey_result <- TukeyHSD(aov_model)
      print(tukey_result)
    }}}

# Perform ANOVA with post-hoc test on each category
AnovaTukey("Area")
AnovaTukey("Assemblage")
AnovaTukey("Type")
AnovaTukey("Preparation")

# Perform MANOVA for PC3 with Type and area/region as interaction
manova_pc3 <- aov(Comp3 ~ Type * Area, data = coord_allpatch)
summary(manova_pc3)


## -- PREF SCARS --
# Define list of PCs and categorical variables
pcs <- paste0("Comp", 1:4)
variables <- c("Area", "Type", "Assemblage", "Preparation")

## Function to perform ANOVA and Tukey post-hoc tests on each PC (1-4)
AnovaTukey <- function(var) {
  for (pc in pcs) {
    formula <- as.formula(paste(pc, "~", var))
    aov_model <- aov(formula, data = coord_allprefs)
    anova_result <- Anova(aov_model, type = "III")
    print(anova_result)
    if (length(unique(coord_allprefs[[var]])) > 2) {
      tukey_result <- TukeyHSD(aov_model)
      print(tukey_result)
    }}}

# Perform ANOVA with post-hoc test on each category
AnovaTukey("Area")
AnovaTukey("Assemblage")
AnovaTukey("Type")
AnovaTukey("Preparation")

# Perform MANOVA for PC1 with Type and area/region as interaction
manova_pc1 <- aov(Comp1 ~ Type * Area, data = coord_allprefs)
summary(manova_pc1)


## -- PRODUCTS --
# Define list of PCs and categorical variables
pcs <- paste0("Comp", 1:6)
variables <- c("Assemblage", "Flake_scars")

## Function to perform ANOVA and Tukey post-hoc tests on each PC (1-4)
AnovaTukey <- function(var) {
  for (pc in pcs) {
    formula <- as.formula(paste(pc, "~", var))
    aov_model <- aov(formula, data = coord_prods)
    anova_result <- Anova(aov_model, type = "III")
    print(anova_result)
    if (length(unique(coord_prods[[var]])) > 2) {
      tukey_result <- TukeyHSD(aov_model)
      print(tukey_result)
    }}}

# Perform ANOVA with post-hoc test on each category
AnovaTukey("Assemblage")
AnovaTukey("Flake_scars")


## -- PREF SCAR AND PRODUCT VENTRAL --
## Function to perform ANOVA and Tukey post-hoc tests on each PC for prefs and prod scars
# Define list of PCs and categorical variables
pcs <- paste0("Comp", 1:4)
variables <- c("Artefact", "Scars", "Assemblage")

AnovaTukey <- function(var) {
  for (pc in pcs) {
    formula <- as.formula(paste(pc, "~", var))
    aov_model <- aov(formula, data = coord_prod_scar)
    anova_result <- Anova(aov_model, type = "III")
    print(anova_result)
    if (length(unique(coord_prod_scar[[var]])) > 2) {
      tukey_result <- TukeyHSD(aov_model)
      print(tukey_result)
    }}}

# Perform ANOVA with post-hoc test on each category
AnovaTukey("Artefact")
AnovaTukey("Scars")
AnovaTukey("Assemblage")

# Perform MANOVA for PC1 and PC3 with assemblage and artefact as interaction
manova_pc1 <- aov(Comp1 ~ Assemblage * Artefact, data = coord_prod_scar)
summary(manova_pc1)

manova_pc3 <- aov(Comp3 ~ Assemblage * Artefact, data = coord_prod_scar)
summary(manova_pc3)

###---------------------------------------------------------------------------

## STEP 8.9 Conduct two-block partial least squares analysis to compare covariation in core outline and preferential scar shape (two.b.pls in geomorph)
# Ensure the core names and core outline/pref specimens match
outlier_spec <- c("CoreME78.627Ap_LR") # this was isolated as an outlier
cores_filtered <- cores[!cores$File_name %in% outlier_spec, ]
specimen_names <- dimnames(alloutline)[[3]]

# Create groups for separate analysis and extract specimens based on Area
NKspecs <- cores_filtered$File_name[cores_filtered$Area == "NK"]
THspecs <- cores_filtered$File_name[cores_filtered$Area == "TH"]

NKoutline <- alloutline[,,specimen_names %in% NKspecs]
THoutline <- alloutline[,,specimen_names %in% THspecs]
NKprefs <- allprefs[,,specimen_names %in% NKspecs]
THprefs <- allprefs[,,specimen_names %in% THspecs]

# Create groups for separate analysis and extract specimens based on Type
T1specs <- cores_filtered$File_name[cores_filtered$Type == "T1"]
T2specs <- cores_filtered$File_name[cores_filtered$Type == "T2"]
T12specs <- cores_filtered$File_name[cores_filtered$Type == "T1/2"]

T1outline <- alloutline[,,specimen_names %in% T1specs]
T2outline <- alloutline[,,specimen_names %in% T2specs]
T12outline <- alloutline[,,specimen_names %in% T12specs]
T1prefs <- allprefs[,,specimen_names %in% T1specs]
T2prefs <- allprefs[,,specimen_names %in% T2specs]
T12prefs <- allprefs[,,specimen_names %in% T12specs]


# Conduct GPA and 2B PLS on the sets of shapes for all specimens
PLS_all <-two.b.pls(outline_gpa$coords, pref_gpa$coords,iter=999)
summary(PLS_all)

# Calculate the percentage of covariation
PLS_all_cov <- (PLS_all$svd$d[1]^2 / sum(PLS_all$svd$d^2)) * 100
print(PLS_all_cov)

# Extract PLS scores
PLS_all_Xscores <- PLS_all$XScores[, 1]  # Scores for Block 1 (outline)
PLS_all_Yscores <- PLS_all$YScores[, 1]  # Scores for Block 2 (preferential scar shape)

# Conduct linear regression analysis on PLS scores
lmPLS_all <- lm(PLS_all_Xscores ~ PLS_all_Yscores)
summary(lmPLS_all)

# Create a new data frame for PLS scores with specimen names and attributes
PLS_all_XY <- data.frame(Specimen = specimen_names, Xscores = PLS_all_Xscores, Yscores = PLS_all_Yscores)
PLS_all_XY <- merge(PLS_all_XY, cores_filtered, by.x = "Specimen", by.y = "File_name")

# Plot X and Y scores by region and add regression line
ggplot(data = PLS_all_XY, aes(x = PLS_all_Xscores, y = PLS_all_Yscores)) +
  geom_point(aes(color = Area, shape = Area), size = 2) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Area), linetype = "dashed", linewidth = 0.7, fullrange = TRUE) +
  ggtitle("All cores: 2B PLS") +
  labs(x = "Block 1/Core outline", y = "Block 2/Pref. scar") +
  theme_minimal() +
  scale_color_manual(values = colsArea) +
  scale_shape_manual(values = shapesArea) +
  geom_hline(yintercept = -0.4, linetype = "solid", color = "darkgrey") + # Horizontal axis line
  geom_vline(xintercept = -0.2, linetype = "solid", color = "darkgrey") + # Vertical axis line
  theme(legend.position = "bottom")

# Plot X and Y scores by core type and add regression line
colsType <- c("T1" = "#FC8D62", "T2" = "#66C2A5", "T1/2" = "#8DA0CB")

ggplot(data = PLS_all_XY, aes(x = PLS_all_Xscores, y = PLS_all_Yscores)) +
  geom_point(aes(color = Type, shape = Area), size = 2) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Type), linetype = "dashed", linewidth = 0.7, fullrange = TRUE) +
  ggtitle("All cores: 2B PLS") +
  labs(x = "Block 1/Core outline", y = "Block 2/Pref. scar") +
  theme_minimal() +
  scale_color_manual(values = colsType) +
  scale_shape_manual(values = shapesArea) +
  geom_hline(yintercept = -0.4, linetype = "solid", color = "darkgrey") + # Horizontal axis line
  geom_vline(xintercept = -0.2, linetype = "solid", color = "darkgrey") + # Vertical axis line
  theme(legend.position = "bottom")

# Repeat analysis for sub-groups
# NK
gpaNKoutline <- gpagen(NKoutline,PrinAxes=TRUE)
gpaNKprefs <- gpagen(NKprefs,PrinAxes=TRUE)

PLS_NK <-two.b.pls(gpaNKoutline$coords, gpaNKprefs$coords,iter=999)
summary(PLS_NK)
PLS_NK_cov <- (PLS_NK$svd$d[1]^2 / sum(PLS_NK$svd$d^2)) * 100
print(PLS_NK_cov)

PLS_NK_Xscores <- PLS_NK$XScores[, 1]
PLS_NK_Yscores <- PLS_NK$YScores[, 1]

lmPLS_NK <- lm(PLS_NK_Xscores ~ PLS_NK_Yscores)
summary(lmPLS_NK)


# TH
gpaTHoutline <- gpagen(THoutline,PrinAxes=TRUE)
gpaTHprefs <- gpagen(THprefs,PrinAxes=TRUE)

PLS_TH <-two.b.pls(gpaTHoutline$coords, gpaTHprefs$coords,iter=999)
summary(PLS_TH)
PLS_TH_cov <- (PLS_TH$svd$d[1]^2 / sum(PLS_TH$svd$d^2)) * 100
print(PLS_TH_cov)

PLS_TH_Xscores <- PLS_TH$XScores[, 1]
PLS_TH_Yscores <- PLS_TH$YScores[, 1]

lmPLS_TH <- lm(PLS_TH_Xscores ~ PLS_TH_Yscores)
summary(lmPLS_TH)


# Type 1
gpaT1outline <- gpagen(T1outline,PrinAxes=TRUE)
gpaT1prefs <- gpagen(T1prefs,PrinAxes=TRUE)

PLS_T1 <-two.b.pls(gpaT1outline$coords, gpaT1prefs$coords,iter=999)
summary(PLS_T1)
PLS_T1_cov <- (PLS_T1$svd$d[1]^2 / sum(PLS_T1$svd$d^2)) * 100
print(PLS_T1_cov)

PLS_T1_Xscores <- PLS_T1$XScores[, 1]
PLS_T1_Yscores <- PLS_T1$YScores[, 1]

lmPLS_T1 <- lm(PLS_T1_Xscores ~ PLS_T1_Yscores)
summary(lmPLS_T1)


# Type 2
gpaT2outline <- gpagen(T2outline,PrinAxes=TRUE)
gpaT2prefs <- gpagen(T2prefs,PrinAxes=TRUE)

PLS_T2 <-two.b.pls(gpaT2outline$coords, gpaT2prefs$coords,iter=999)
summary(PLS_T2)
PLS_T2_cov <- (PLS_T2$svd$d[1]^2 / sum(PLS_T2$svd$d^2)) * 100
print(PLS_T2_cov)

PLS_T2_Xscores <- PLS_T2$XScores[, 1]
PLS_T2_Yscores <- PLS_T2$YScores[, 1]

lmPLS_T2 <- lm(PLS_T2_Xscores ~ PLS_T2_Yscores)
summary(lmPLS_T2)


# Type 1/2
gpaT12outline <- gpagen(T12outline,PrinAxes=TRUE)
gpaT12prefs <- gpagen(T12prefs,PrinAxes=TRUE)

PLS_T12 <-two.b.pls(gpaT12outline$coords, gpaT12prefs$coords,iter=999)
summary(PLS_T12)
PLS_T12_cov <- (PLS_T12$svd$d[1]^2 / sum(PLS_T12$svd$d^2)) * 100
print(PLS_T12_cov)

PLS_T12_Xscores <- PLS_T12$XScores[, 1]
PLS_T12_Yscores <- PLS_T12$YScores[, 1]

lmPLS_T12 <- lm(PLS_T12_Xscores ~ PLS_T12_Yscores)
summary(lmPLS_T12)

###---------------------------------------------------------------------------

## STEP 8.10 Assess shape standardisation using Coefficient of Variation on Procrustes distances
## -- CORE SURFACES --
# Extract gpa coords and calculate mean shape and Procrustes distances from the mean shape
patchcoords <- patched_gpa$coords
patchmeansh <- mshape(patchcoords)
patchproc_dist <- apply(patchcoords, 3, function(x) sqrt(sum((x - patchmeansh)^2)))

# Associate specimen names with Procrustes distance data
specimen_names <- dimnames(patchcoords)[[3]]
patchproc_dist_tab <- data.frame(Specimen = specimen_names, Distance = patchproc_dist)

patchproc_dist_tab <- patchproc_dist_tab %>%
  left_join(cores_filtered, by = c("Specimen" = "File_name"))

# Calculate CV (mean/sd * 100) by area/region and assemblage
cv_patch_Area <- patchproc_dist_tab %>%
  #mutate(Area = ifelse(is.na(Area), "NK", Area)) %>%
  group_by(Area) %>%
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV = sd(Distance) / mean(Distance) * 100)
print(cv_patch_Area)

cv_patch_Assemblage <- patchproc_dist_tab %>%
  #mutate(Assemblage = ifelse(is.na(Assemblage), "NK1_M", Assemblage)) %>%
  group_by(Assemblage) %>%
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV = sd(Distance) / mean(Distance) * 100)
print(cv_patch_Assemblage)


## -- PREF SCARS --
# Extract gpa coords and calculate mean shape and Procrustes distances from the mean shape
prefcoords <- pref_gpa$coords
prefmeansh <- mshape(prefcoords)
prefproc_dist <- apply(prefcoords, 3, function(x) sqrt(sum((x - prefmeansh)^2)))

# Associate specimen names with Procrustes distance data
specimen_names <- dimnames(prefcoords)[[3]]
prefproc_dist_tab <- data.frame(Specimen = specimen_names, Distance = prefproc_dist)

prefproc_dist_tab <- prefproc_dist_tab %>%
  left_join(cores_filtered, by = c("Specimen" = "File_name"))

# Calculate CV (mean/sd * 100) by area/region and assemblage
cv_pref_Area <- prefproc_dist_tab %>%
  group_by(Area) %>%
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV_pd = sd(Distance) / mean(Distance) * 100)  # CV as percentage
print(cv_pref_Area)

cv_pref_Assemblage <- prefproc_dist_tab %>%
  group_by(Assemblage) %>%
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV_pd = sd(Distance) / mean(Distance) * 100) # CV as percentage
print(cv_pref_Assemblage)


## -- NK PRODUCT OUTLINES --
# Extract gpa coords and calculate mean shape and Procrustes distances from the mean shape
prodvcoords <- prodv_gpa$coords
prodvmeansh <- mshape(prodvcoords)
prodvproc_dist <- apply(prodvcoords, 3, function(x) sqrt(sum((x - prodvmeansh)^2)))

# Associate specimen names with Procrustes distance data
specimen_names <- dimnames(prodvcoords)[[3]]
prodvproc_dist_tab <- data.frame(Specimen = specimen_names, Distance = prodvproc_dist)

prodvproc_dist_tab <- prodvproc_dist_tab %>%
  left_join(prods_filtered, by = c("Specimen" = "File_name"))

# Calculate CV (mean/sd * 100) by area/region and assemblage
cv_prodv_Assemblage <- prodvproc_dist_tab %>%
  group_by(Assemblage) %>%
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV = sd(Distance) / mean(Distance) * 100)
print(cv_prodv_Assemblage)

cv_prodv_NK <- prodvproc_dist_tab %>%  # overall CV for NK products
  summarise(
    Mean_pd = mean(Distance),
    SD_pd = sd(Distance),
    CV = sd(Distance) / mean(Distance) * 100)
print(cv_prodv_NK)


## EXTENDED METHODS STEP 5A: R SCRIPT FOR PREPARING CORE DATA
### Alternative method without distal landmark

# Install and load required packages
library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(dplyr)
library(tibble)
library(ggplot2)

###---------------------------------------------------------------------------
## STEP 5.1 Extract core LMs from .pts and .ply files
# Define a path to the folder with core data
core_data <- file.path("analysis/data/raw_data/core_data")

# Create list of core names
corelist <- list.files(path = core_data, pattern = ".pts")

# Get names from files
names <- gsub(".pts", "", basename(corelist))

# Create empty array
nd.core_landmarks <- NULL

# Import each file (selecting landmarks 2 to 63) (i.e. not the distal LM) and paste them in a single matrix
for (i in 1:length(corelist)){
  core_path <- file.path(core_data, corelist[i])
  tmp <- as.data.frame(read.table(core_path, skip = 2, header = F)[2:63,2:4])
  tmp <- as.matrix(tmp)
  nd.core_landmarks <- rbind(nd.core_landmarks, tmp)
}

nd.core_landmarks <- mapply(nd.core_landmarks, FUN=as.numeric)
nd.core_landmarks <- matrix(data=nd.core_landmarks, ncol=3)

# Convert landmark 3D data matrix into array (62 lms, 3 dimensions)
nd.core_landmarks <- geomorph::arrayspecs(nd.core_landmarks, 62, 3)

# View specimen example with imported landmark configuration
core1 <- vcgImport(file.path(core_data, "CoreME78.627Aa_LR.ply"))
shade3d(core1, color = "gray")
spheres3d(nd.core_landmarks[1:2,,34], color = 1, radius = 1)
spheres3d(nd.core_landmarks[3:22,,34], color = 2, radius = 0.7)
spheres3d(nd.core_landmarks[23:62,,34], color = 3, radius = 0.7)
#rgl.snapshot('nd.core_landmarks_example.png', fmt = 'png')
#close3d()

###---------------------------------------------------------------------------

## STEP 5.2 Resample core outline semilandmarks
# Create a new data array for resampled semilandmarks (38 slms, 3 dimensions, 166 specimens)
nd.core_lms <- array(0, dim=c(38,3,166))

# Resample semilandmarks on curves between fixed landmarks and redistribute to be equally-spaced
for(i in 1:dim(nd.core_landmarks)[3]) {
  startA <- nd.core_landmarks[2,,i]
  startB <- nd.core_landmarks[1,,i]
  curveA <- digit.curves(startA, rbind(nd.core_landmarks[3:22,,i], startB), nPoints = 6, closed = F)
  curveB <- digit.curves(startB, rbind(nd.core_landmarks[23:62,,i], startA), nPoints = 30, closed = F)
  comb <- rbind(startA, startB, curveA[2:(nrow(curveA)-1),], curveB[2:(nrow(curveB)-1),])
  nd.core_lms[,,i] <- comb
}

# Save the new semilandmark array for ease of future use
save(nd.core_lms, file = "nd.core_lms.RData")

# To view an example core mesh (.ply) with resampled slms (identify the correct specimen number (#34) by viewing 'corelist')
shade3d(core1, color = "gray")
spheres3d(nd.core_lms[1:2,,34], color = 1, radius = 1)
spheres3d(nd.core_lms[3:8,,34], color = 2, radius = 0.7)
spheres3d(nd.core_lms[9:38,,34], color = 3, radius = 0.7)
#rgl.snapshot('nd.core_lms_example.png', fmt = 'png')
#close3d()

## STEP 5.3 Relax semilandmarks on curve outlines
# Define fixed LMs (that do not slide) and semilandmarks on outline curves to slide
fix <- c(1:2)
outlineA <- c(3:8)
outlineB <- c(9:38)

outlines <- list(outlineA, outlineB)

dimnames(nd.core_lms)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curves <- slider3d(nd.core_lms, SMvector = fix, deselect = TRUE, outlines = outlines,
                         sur.path = "core_data", sur.type = "ply", iterations = 3)

# Visualise slid slms from previous to new positions
deformGrid3d(slide_curves$dataslide[,,34],nd.core_lms[,,34],ngrid = 0)

# View slid slms on mesh model
shade3d(core1, color = "gray")
spheres3d(slide_curves$dataslide[fix,,34],col=1,radius=1)
spheres3d(slide_curves$dataslide[,,34],col=3,radius=0.7)
#rgl.snapshot('nd.core_lms_slid.png', fmt = 'png')
#close3d()

# Rename and save outline slms
nd.alloutline <- slide_curves$dataslide
save(nd.alloutline, file = "nd.alloutline.RData")

###--Step 5.4 is unaltered----------------------------------------------------

## STEP 5.4 Extract preferential scar outline on all specimens
# Create empty array
pref_landmarks <- NULL

# Import each file and paste them in a single matrix (this extracts only landmarks 71 to 120 from the .pts file)
for (i in 1:length(corelist)){
  core_path <- file.path(core_data, corelist[i])
  tmp <- as.data.frame(read.table(core_path, skip = 2, header = F)[71:120,2:4])
  tmp <- as.matrix(tmp)
  pref_landmarks <- rbind(pref_landmarks, tmp)
}

pref_landmarks <- mapply(pref_landmarks, FUN=as.numeric)
pref_landmarks <- matrix(data=pref_landmarks, ncol=3)

# Convert landmark 3D data matrix into array (50 lms, 3 dimensions)
pref_landmarks <- geomorph::arrayspecs(pref_landmarks, 50, 3)

# View specimen example with imported landmark configuration
shade3d(core1, color = "gray")
spheres3d(pref_landmarks[1,,34], color = 1, radius = 1)
spheres3d(pref_landmarks[40,,34], color = 1, radius = 1)
spheres3d(pref_landmarks[2:39,,34], color = 3, radius = 0.7)
spheres3d(pref_landmarks[41:50,,34], color = 4, radius = 0.7)
#rgl.snapshot('pref_landmarks_example.png', fmt = 'png')

# Create an empty array for resampled semilandmarks (37 slms, 3 dimensions, 166 specimens)
pref_lms <- array(0, dim=c(37,3,166))

# Redistribute equally-spaced semilandmarks on all specimens along two curves (5 platform and 30 edge points)
for(i in 1:dim(pref_landmarks)[3]) {
  startA <- pref_landmarks[1,,i]
  startB <- pref_landmarks[40,,i]
  curveA <- digit.curves(startA,rbind(pref_landmarks[2:39,,i], startB), nPoints = 30, closed = F)
  curveB <- digit.curves(startB, rbind(pref_landmarks[41:50,,i], startA), nPoints = 5, closed = F)
  comb <- rbind(curveA[2:nrow(curveA),], curveB[2:nrow(curveB),])
  pref_lms[,,i] <- comb
}

shade3d(core1, color = "gray")
spheres3d(pref_lms[37,,34], color = 1, radius = 1)
spheres3d(pref_lms[31,,34], color = 1, radius = 1)
spheres3d(pref_lms[1:30,,34], color = 3, radius = 0.7)
spheres3d(pref_lms[32:36,,34], color = 4, radius = 0.7)
#rgl.snapshot('pref_lms_example.png', fmt = 'png')

# Save for future use
save(pref_lms, file = "pref_lms.RData")

# Fix platform landmarks and relax edge semilandmarks on pref scar outline
fixP <- c(31:37)
outlinesP <- c(1:30)

dimnames(pref_lms)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curvesP <- slider3d(pref_lms, SMvector = fixP, deselect = TRUE, outlines = outlinesP,
                          sur.path = "core_data", sur.type = "ply", iterations = 3)

# Visualise slid slms from previous to new positions
deformGrid3d(slide_curvesP$dataslide[,,34],pref_lms[,,34],ngrid = 0)

# View slid slms on mesh model
shade3d(core1, color = "gray")
spheres3d(slide_curvesP$dataslide[fixP,,34],col=4,radius=1)
spheres3d(slide_curvesP$dataslide[,,34],col=3,radius=0.7)
#rgl.snapshot('pref_lms_slid.png', fmt = 'png')

# Rename and save
allprefs <- slide_curvesP$dataslide
save(allprefs, file = "allprefs.RData")

###--Resume alterations here---------------------------------------------------

## STEP 5.5 Use semilandmarked surface patches on a single specimen to create a template
# Create atlas using (arbitrarily selected) core ME78.389c (specimen #20 for slid outline slms) as a template
coreT <- vcgImport(file.path(core_data, "CoreME78.389c_LR.ply"))
core_temp <- as.matrix(read.table("CoreME78.389c_LRtemplate.pts", skip = 2, header = F)[,2:4])

nd.core_atlas <- createAtlas(coreT, landmarks = nd.alloutline[,,20], patch = core_temp)
plotAtlas(nd.core_atlas)
#close3d()

###---------------------------------------------------------------------------

## STEP 5.6 Place the patch onto all specimens. This deforms the template over each mesh surface. The inflate parameter can be adjusted if there are placement errors observed in the checking stage
nd.allpatched <- placePatch(atlas = nd.core_atlas, dat.array = nd.alloutline, path = "core_data", inflate = 10, fileext = ".ply")

# Save the patched slm data
save(nd.allpatched, file = "nd.allpatched.RData")

###---------------------------------------------------------------------------

## STEP 5.7 Check patch slms have deformed correctly on each specimen (this is time consuming but essential to check the deformation process). If there are placement errors, try increasing the inflate parameter, or repeat previous steps excluding problem specimens.
checkLM(nd.allpatched, atlas=nd.core_atlas)

## Problems in deformation identified in specimens 53, 59, 67, 82, 88, 104, 128 (minor) and 87, 99, 125, 133, 134 (major)

# Separate the failed specimens
failed <- c(53, 59, 67, 82, 88, 104, 128, 87, 99, 125, 133, 134)

# Separate failed and successful specimens
exclude <- rep(FALSE, 166)
exclude[failed] <- TRUE

outline_failed <- nd.alloutline[,,exclude]
patched_good <- nd.allpatched[,,!exclude]


# Re-run the patch deformation on the failed specimens with a higher inflate parameter
nd.patched_new <- placePatch(atlas = nd.core_atlas, dat.array = outline_failed, path = "core_data", inflate = 10, fileext = ".ply")
checkLM(nd.patched_new, atlas=nd.core_atlas)

## Problems still identified in specimens 2, 3, 5, 7, 9, 11, 12 - repeat the exclusion process
# Separate the failed specimens
failed <- c(2, 3, 5, 7, 9, 11, 12)

# Separate failed and successful specimens
exclude <- rep(FALSE, 12)
exclude[failed] <- TRUE
patched_new2 <- nd.patched_new[,,!exclude]

# Create a new array and combine both sets of successfully patched specimens on the 3rd dimension
nd.allpatched <- array(NA, dim = c(dim(patched_good)[1], dim(patched_new2)[2], 159))
nd.allpatched <- abind(patched_good, patched_new2, along = 3)

# Ensure the dimnames are preserved
dimnames(nd.allpatched)[[3]] <- c(dimnames(patched_good)[[3]], dimnames(patched_new2)[[3]])


## STEP 5.8 Check for errors in alignment (e.g. inversions) by conducting GPA and inspecting landmark configurations.
#These errors will also be apparent as outliers when plotted on a PCA.
nd.patched_gpa <- gpagen(nd.allpatched,PrinAxes=TRUE)

# Extract the coordinate data from the gpa
coords <- nd.patched_gpa$coords

# Get the number of specimens
num_specimens <- dim(coords)[3]

# Loop through each specimen and create a separate plot, then loop through each landmark to plot configurations in different colours. These are saved individually.
for (specimen_idx in 1:num_specimens) {
  specimen_coords <- coords[, , specimen_idx]
  for (landmark_idx in 1:438) {
    landmark <- specimen_coords[landmark_idx, , drop = FALSE]
    x <- as.numeric(landmark[, 1]) # Extract x,y,z coordinates
    y <- as.numeric(landmark[, 2])
    z <- as.numeric(landmark[, 3])
    if (landmark_idx %in% c(1:2)) {
      col <- "red"  # Fixed LM1-2 in red
    } else if (landmark_idx %in% c(3:38)) {
      col <- "blue"  # Outline LMs3-38 in blue
    } else if (landmark_idx %in% c(39:238)) {
      col <- "green"  # Upper surface patch in green
    } else {
      col <- "yellow"  # Lower surface patch in yellow
    }
    points3d(x, y, z, col = col, size = 5)
  }
  view3d(theta = 180, phi = 0, zoom = 1, fov = 0, interactive = TRUE)
  rgl.snapshot(file = paste("specimen_", specimen_idx, ".png", sep = ""), fmt = "png")
  close3d()
}

## 22 specimens show inversions: 10, 14, 31, 35, 37, 38, 62, 67, 69, 73, 76, 78, 80, 97, 104, 111, 143, 145, 155, 156, 159

# Identify and exclude the failed specimen numbers
failed <- c(10, 14, 31, 35, 37, 38, 62, 67, 69, 73, 76, 78, 80, 97, 104, 111, 143, 145, 155, 156, 159)
nd.filtered <- nd.allpatched[, , -failed]

# Repeat GPA
nd.patched_gpa <- gpagen(nd.filtered,PrinAxes=TRUE)


### Modified core data are now ready for GM analyses

## R data files can be loaded as:
load("nd.allpatched.RData")
load("nd.alloutline.RData")


################################################################

## STEP 8.3 Conduct PCAs (gm.prcomp in geomorph)

# Check for outliers:
# Calculate Procrustes distances to the mean shape
meansh <- mshape(nd.patched_gpa$coords)
proc_dist <- apply(nd.patched_gpa$coords, 3, function(x) sqrt(sum((x - meansh)^2)))

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

# Exclude outlier spec (CoreME78.373.37) from the (filtered) patched core dataset and rerun GPA on new dataset
outlier_specs <- 56
all_specs <- 1:138
keep <- setdiff(all_specs, outlier_specs)
nd.filtered <- nd.filtered[, , keep]
nd.patched_gpa <- gpagen(nd.filtered,PrinAxes=TRUE) # check this has the right no of specimens


# Conduct PCA and extract PC scores to new dataframe
PCA_allpatch.nd <- gm.prcomp(nd.patched_gpa$coords)
PCscores_patch.nd <- PCA_allpatch.nd$x # x represents the PC scores
coord_allpatch.nd <- as.data.frame(PCscores_patch.nd)

# Import attribute files
cores <- read.csv("File_list_cores.csv")
cores <- dplyr::arrange(cores, File_name)

# Match the coord to the core attribute data only for the specimens remaining in the sample
coord_allpatch.nd$File_name <- rownames(coord_allpatch.nd)
coord_allpatch.nd <- merge(coord_allpatch.nd, cores, by = "File_name")


# Generate variance tables that show the percentage of variation captured by each component
eigenvalues <- PCA_allpatch.nd$d # d represents the eigenvalues
pc_numbers <- 1:length(eigenvalues)
perc_variance <- (eigenvalues / sum(eigenvalues)) * 100
cum_perc <- cumsum(perc_variance)
perc_table_patch.nd <- data.frame(PC = 1:length(eigenvalues), Eigenvalue = eigenvalues, PercVariance = perc_variance, CumulativePerc = cum_perc)


############## PLOT PLCA

## STEP 8.4 Plot PCAs. Repeat for products and showing other variables (e.g. core preparation) as aesthetics

## Generate PCA plots (here with aesthetics set for area/region)
# Set colour/shape aesthetics for the NK and TH samples for the plots
colsArea <- c("NK" = "royalblue", "TH" = "indianred1")
shapesArea <- c("NK" = 16, "TH" = 17)

#PC1 and 2
coord_allpatch.nd %>%
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

#PC3 and 4
#CoreME78.737.37_LR (Spec 56) exhibits anomalous values on PC3
coord_allpatch.nd %>%
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
coord_allpatch.nd %>%
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

################################

## STEP 8.5 Additional visualisation steps to represent shape variation of core surfaces
# Find mean, minimum and maximum shapes (mshape and PlotReftoTarget visualisation in geomorph)

## -- CORE SURFACES --
# Find the theoretical mean shape
msh_patch <- mshape(nd.patched_gpa$coords)

# Create an outline by linking consecutive slms of the shape
#define.links(msh_patch, ptsize = 3, links = NULL) # can define links between slms manually or input as a .csv
core_links.nd <- read.csv("Core_links.nd.csv")

# Identify the specimen that most closely matches the theoretical mean shape
findMeanSpec(nd.patched_gpa$coords)
#CoreME78.778d_LR
#64

# Define a function to generate and save min./max. shape snapshots for a given PC
GeneratePCShapes <- function(pc_num) {
  # Mean shape (grey) and minimum shape (blue)
  plotRefToTarget(msh_patch, PCA_allpatch.nd$shapes[[paste0("shapes.comp", pc_num)]]$min,
                  method = "points", links = core_links.nd,
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
  plotRefToTarget(msh_patch, PCA_allpatch.nd$shapes[[paste0("shapes.comp", pc_num)]]$max,
                  method = "points", links = core_links.nd,
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
coord_allpatch.nd <- data.frame(
  Specimen = coord_allpatch.nd$File_name,
  coord_allpatch.nd,
  stringsAsFactors = FALSE)

# Extract specimen names that correspond to min. and max. value for first 6 PCs
for (i in 1:6) {
  pc_name <- paste0("Comp", i)
  min_specimen <- coord_allpatch.nd$Specimen[which.min(coord_allpatch.nd[[pc_name]])]
  max_specimen <- coord_allpatch.nd$Specimen[which.max(coord_allpatch.nd[[pc_name]])]
  Specs <- rbind(Specs, data.frame(PC = pc_name, MinSp = min_specimen, MaxSp = max_specimen, stringsAsFactors = FALSE))
}
print(Specs)


## Continue with outher STEP 8 GM Analyses as needed.

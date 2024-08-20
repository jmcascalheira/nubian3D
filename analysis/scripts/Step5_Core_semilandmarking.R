
## STEP 5. R SCRIPT FOR PREPARING CORE DATA

# Load required packages
library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(dplyr)

###------------------------------------------------------------------------
## STEP 5.1 Extract core LMs from .pts and .ply files
# Define a path to the folder with core data
core_data <- file.path("analysis/data/raw_data/core_data")

# Create list of core names
corelist <- list.files(path = core_data, pattern = ".pts")

# Get names from files
names <- gsub(".pts", "", basename(corelist))

# Create empty array
core_landmarks <- NULL

# Import each file (selecting landmarks 1 to 63) and paste them in a single matrix
for (i in 1:length(corelist)){
  core_path <- file.path(core_data, corelist[i])
  tmp <- as.data.frame(read.table(core_path, skip = 2, header = F)[1:63,2:4])
  tmp <- as.matrix(tmp)
  core_landmarks <- rbind(core_landmarks, tmp)
}

core_landmarks <- mapply(core_landmarks, FUN=as.numeric)
core_landmarks <- matrix(data=core_landmarks, ncol=3)

# Convert landmark 3D data matrix into array (63 lms, 3 dimensions)
core_landmarks <- geomorph::arrayspecs(core_landmarks, 63, 3)

# View specimen example with imported landmark configuration
core1 <- vcgImport(file.path(core_data, "CoreME78.627Aa_LR.ply"))
shade3d(core1, color = "gray")
spheres3d(core_landmarks[1:3,,34], color = 1, radius = 1)
spheres3d(core_landmarks[4:23,,34], color = 2, radius = 0.7)
spheres3d(core_landmarks[24:43,,34], color = 3, radius = 0.7)
spheres3d(core_landmarks[44:63,,34], color = 4, radius = 0.7)
#rgl.snapshot('core_landmarks_example.png', fmt = 'png')
#close3d()

###------------------------------------------------------------------------
## STEP 5.2 Resample core outline semilandmarks
# Create a new data array for resampled semilandmarks (39 slms, 3 dimensions, 166 specimens)
core_lms <- array(0, dim=c(39,3,166))

# Resample semilandmarks on curves between fixed landmarks and redistribute to be equally-spaced
for(i in 1:dim(core_landmarks)[3]) {
  startA <- core_landmarks[3,,i]
  startB <- core_landmarks[2,,i]
  startC <- core_landmarks[1,,i]
  curveA <- digit.curves(startA, rbind(core_landmarks[4:23,,i], startB), nPoints = 6, closed = F)
  curveB <- digit.curves(startB, rbind(core_landmarks[24:43,,i], startC), nPoints = 15, closed = F)
  curveC <- digit.curves(startC, rbind(core_landmarks[44:63,,i], startA), nPoints = 15, closed = F)
  comb <- rbind(startA, startB, startC, curveA[2:(nrow(curveA)-1),], curveB[2:(nrow(curveB)-1),], curveC[2:(nrow(curveC)-1),])
  core_lms[,,i] <- comb
}

# Save the new semilandmark array for ease of future use
save(core_lms, file = "analysis/data/derived_data/core_lms.RData")

# To view an example core mesh (.ply) with resampled slms (identify the correct specimen number (#34) by viewing 'corelist')
shade3d(core1, color = "gray")
spheres3d(core_lms[1:3,,34], color = 1, radius = 1)
spheres3d(core_lms[4:9,,34], color = 2, radius = 0.7)
spheres3d(core_lms[10:24,,34], color = 3, radius = 0.7)
spheres3d(core_lms[25:39,,34], color = 4, radius = 0.7)
#rgl.snapshot('core_lms_example.png', fmt = 'png')
#close3d()

###---------------------------------------------------------------------------

## STEP 5.3 Relax semilandmarks on curve outlines
# Define fixed LMs (that do not slide) and semilandmarks on outline curves to slide
fix <- c(1:3)
outlineA <- c(4:9)
outlineB <- c(10:24)
outlineC <- c(25:39)

outlines <- list(outlineA, outlineB, outlineC)

dimnames(core_lms)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curves <- slider3d(core_lms, SMvector = fix, deselect = TRUE, outlines = outlines,
                         sur.path = "core_data", sur.type = "ply", iterations = 3)

# Visualise slid slms from previous to new positions
deformGrid3d(slide_curves$dataslide[,,34],core_lms[,,34],ngrid = 0)
#close3d()

# View slid slms on mesh model
shade3d(core1, color = "gray")
spheres3d(slide_curves$dataslide[fix,,34],col=1,radius=1)
spheres3d(slide_curves$dataslide[,,34],col=3,radius=0.7)
#rgl.snapshot('core_lms_slid.png', fmt = 'png')
#close3d()

# Rename and save outline slms
alloutline <- slide_curves$dataslide
save(alloutline, file = "analysis/data/derived_data/alloutline.RData")

###---------------------------------------------------------------------------

## STEP 5.4 Extract preferential scar outline on all specimens
# Create a new data array for pref scar landmarks
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
#close3d()

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

# View resampled example specimen
shade3d(core1, color = "gray")
spheres3d(pref_lms[37,,34], color = 1, radius = 1)
spheres3d(pref_lms[31,,34], color = 1, radius = 1)
spheres3d(pref_lms[1:30,,34], color = 3, radius = 0.7)
spheres3d(pref_lms[32:36,,34], color = 4, radius = 0.7)
#rgl.snapshot('pref_lms_example.png', fmt = 'png')
#close3d()

# Save for future use
save(pref_lms, file = "analysis/data/derived_data/pref_lms.RData")

# Fix platform landmarks and relax edge semilandmarks on pref scar outline
fixP <- c(31:37)
outlinesP <- c(1:30)

dimnames(pref_lms)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curvesP <- slider3d(pref_lms, SMvector = fixP, deselect = TRUE, outlines = outlinesP,
                          sur.path = "core_data", sur.type = "ply", iterations = 3)

# Visualise slid slms from previous to new positions
deformGrid3d(slide_curvesP$dataslide[,,34],pref_lms[,,34],ngrid = 0)
#close3d()

# View slid slms on mesh model
shade3d(core1, color = "gray")
spheres3d(slide_curvesP$dataslide[fixP,,34],col=4,radius=1)
spheres3d(slide_curvesP$dataslide[,,34],col=3,radius=0.7)
#rgl.snapshot('pref_lms_slid.png', fmt = 'png')
#close3d()

# Rename and save
allprefs <- slide_curvesP$dataslide
save(allprefs, file = "analysis/data/derived_data/allprefs.RData")

###---------------------------------------------------------------------------

## STEP 5.5 Use semilandmarked surface patches on a single specimen to create a template
# Create atlas using (arbitrarily selected) core ME78.389c (specimen #20 for slid outline slms) as a template
coreT <- vcgImport(file.path(core_data, "CoreME78.389c_LR.ply"))
core_temp <- as.matrix(read.table("CoreME78.389c_LRtemplate.pts", skip = 2, header = F)[,2:4])

core_atlas <- createAtlas(coreT, landmarks = alloutline[,,20], patch = core_temp)
plotAtlas(core_atlas)
#close3d()

###---------------------------------------------------------------------------

## STEP 5.6 Place the patch onto all specimens
# Deform the template over each mesh surface. The inflate parameter can be adjusted if there are placement errors observed in the checking stage
allpatched <- placePatch(atlas = core_atlas, dat.array = alloutline, path = "core_data", inflate = 7, fileext = ".ply")

###---------------------------------------------------------------------------

## STEP 5.7 Check patch slms have deformed correctly on each specimen (this is time consuming but essential to check the deformation process). If there are placement errors, try increasing the inflate parameter, or repeat previous steps excluding problem specimens.
checkLM(allpatched, atlas=core_atlas)

# Save the patched slm data
save(allpatched, file = "analysis/data/derived_data/allpatched.RData")

###---------------------------------------------------------------------------

## STEP 5.8 Create separate datasets for NK and TH samples and append attributes to specimens
# Extract specimens 1-115 into "patchedNK"
patchedNK <- allpatched[, , 1:115]
save(patchedNK, file = "analysis/data/derived_data/patchedNK.RData")

# Extract specimens 116-166 into "patchedTH"
patchedTH <- allpatched[, , 116:166]
save(patchedTH, file = "analysis/data/derived_data/patchedTH.RData")

# Extract specimens 1-115 into "prefNK"
prefNK <- allprefs[, , 1:115]
save(prefNK, file = "analysis/data/derived_data/prefNK.RData")

# Extract specimens 116-166 into "prefTH"
prefTH <- allprefs[, , 116:166]
save(prefTH, file = "analysis/data/derived_data/prefTH.RData")


# Import attributes for specimens
cores <- read.csv("File_list_cores.csv")
cores <- dplyr::arrange(cores, File_name)
Assemblage <- factor(cores$Assemblage)
Area <- factor(cores$Area)
Scars <- factor(cores$Preparation)
Type <- factor(cores$Type)

###---------------------------------------------------------------------------

## Core data are now ready for GM analyses

## R data files can be loaded as:
#core surfaces
load("allpatched.RData")
load("patchedNK.RData")
load("patchedTH.RData")

#core outlines only
load("alloutline.RData")

#pref scar outlines only
load("allprefs.RData")
load("prefNK.RData")
load("prefTH.RData")


## STEP 6. R SCRIPT FOR PREPARING PRODUCT (BLANK/FLAKE) DATA

# Ensure all required libraries are loaded
library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(dplyr)

###---------------------------------------------------------------------------

## STEP 6.1 Extract product LMs from .pts and .ply files
# Define a path to the folder with core data
Product_data <- file.path("analysis/data/raw_data/product_data")

# Create list of files with .pts extension
prodlist <- list.files(path = Product_data, pattern = ".pts")

# Get names from files
names <- gsub(".pts", "", basename(prodlist))

# Create empty array
prod_landmarks <- NULL

# Import each file and paste them in a single matrix. By placing points 4-63 in descending order, this step reverses the direction of the landmarks
for (i in 1:length(prodlist)){
  prod_path <- file.path(Product_data, prodlist[i])
  tmp <- as.data.frame(read.table(prod_path, skip = 2, header = F)[,2:4])
  fix <- tmp[1:3,]
  tmp <- arrange(tmp, desc(row_number()))
  tmp <- rbind(fix,tmp)
  tmp <- tmp[1:63,]
  tmp <- as.matrix(tmp)
  prod_landmarks <- rbind(prod_landmarks, tmp)
}

prod_landmarks <- mapply(prod_landmarks, FUN=as.numeric)
prod_landmarks <- matrix(data=prod_landmarks, ncol=3)

# Convert landmark 3D data matrix into array (63 lms, 3 dimensions)
prod_landmarks <- geomorph::arrayspecs(prod_landmarks, 63, 3)

# To check the new landmark configuration on an example mesh (identify the correct specimen number (#3) by viewing 'prodlist')
product1 <- vcgImport(file.path(Product_data, "PointME78.1017.4_LR.ply"))
shade3d(product1, color = "gray")
spheres3d(prod_landmarks[1,,3], color = 1, radius = 1)
spheres3d(prod_landmarks[2:3,,3], color = "darkgrey", radius = 1)
spheres3d(prod_landmarks[4:43,,3], color = 2, radius = 0.7)
spheres3d(prod_landmarks[44:48,,3], color = 3, radius = 0.7)
spheres3d(prod_landmarks[49:58,,3], color = 4, radius = 0.7)
spheres3d(prod_landmarks[59:63,,3], color = 5, radius = 0.7)
#rgl.snapshot('analysis/figures/prod_landmarks_example.png', fmt = 'png')
close3d()

###---------------------------------------------------------------------------

## STEP 6.2 Resample product outline semilandmarks
# Create a new data array for resampled semilandmarks (47 slms, 3 dimensions, 178 specimens)
prod_lms <- array(0, dim=c(47,3,178))

# Resample semilandmarks on curves between fixed landmarks and redistribute to be equally-spaced
for(i in 1:dim(prod_landmarks)[3]) {
  pp <- prod_landmarks[1,,i]
  ptA <- prod_landmarks[3,,i]
  ptB <- prod_landmarks[2,,i]
  curveA <- digit.curves(ptB, rbind(prod_landmarks[59:63,,i], pp), nPoints = 3, closed = F)
  curveB <- digit.curves(pp, rbind(prod_landmarks[44:48,,i], ptA), nPoints = 3, closed = F)
  curveC <- digit.curves(ptA, rbind(prod_landmarks[49:58,,i], ptB), nPoints = 8, closed = F)
  curveD <- digit.curves(ptA,rbind(prod_landmarks[4:43,,i],ptB), nPoints = 30, closed = F)
  comb <- rbind(pp, ptB, ptA, curveA[2:(nrow(curveA)-1),], curveB[2:(nrow(curveB)-1),], curveC[2:(nrow(curveC)-1),], curveD[2:(nrow(curveD)-1),])
  prod_lms[,,i] <- comb
}

# Save the new semilandmark array for ease of future use
save(prod_lms, file = "analysis/data/derived_data/prod_lms.RData")

# To view the resampled slm configuration on an example mesh
shade3d(product1, color = "gray")
spheres3d(prod_lms[1,,3], color = 1, radius = 1) #fixed lms
spheres3d(prod_lms[2:3,,3], color = "darkgrey", radius = 1) #fixed lms
spheres3d(prod_lms[4:6,,3], color = 5, radius = 0.7) #ventr platform slms
spheres3d(prod_lms[7:9,,3], color = 3, radius = 0.7) #ventr platform slms
spheres3d(prod_lms[10:17,,3], color = 4, radius = 0.7) #dors platform slms
spheres3d(prod_lms[18:47,,3], color = 2, radius = 0.7) #outline slms
#rgl.snapshot('analysis/figures/prod_lm_example.png', fmt = 'png')
close3d()

###---------------------------------------------------------------------------

## STEP 6.3 Relax semilandmarks on curve outlines
# Define fixed LMs (that do not slide) and semilandmarks on outline curves to slide
fix <- c(1:3)
outlineA <- c(4:6)
outlineB <- c(7:9)
outlineC <- c(10:17)
outlineD <- c(18:47)

outlines <- list(outlineA, outlineB, outlineC, outlineD)

dimnames(prod_lms)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curves <- slider3d(prod_lms, SMvector = fix, deselect = TRUE, outlines = outlines,
                            sur.path = "analysis/data/raw_data/product_data", sur.type = "ply", iterations = 3)

# Visualise slid slms from previous to new positions
deformGrid3d(slide_curves$dataslide[,,3],prod_lms[,,3],ngrid = 0)
close3d()

# View slid slms on mesh model
shade3d(product1, color = "gray")
spheres3d(slide_curves$dataslide[fix,,3], color = 1, radius=1)
spheres3d(slide_curves$dataslide[,,3], color = 2, radius=0.7)
#close3d()

# Rename and save
prod_outline <- slide_curves$dataslide # Product outline data are not used for analysis in this paper but is used in subsequent slm processing steps
save(prod_outline, file = "analysis/data/derived_data/prod_outline.RData")

###---------------------------------------------------------------------------

## STEP 6.4 Use semilandmarked surface patches on a single specimen to create a template
# Create atlas using (arbitrarily selected) product ME78.1017.4 (specimen #3 for slid outline slms) as a template
prod_temp <- as.matrix(read.table("analysis/data/derived_data/PointME78.1017.4_LRtemplate.pts", skip = 2, header = F)[,2:4])
outl_temp <- prod_outline[,,3]

prod_atlas <- createAtlas(product1, landmarks = outl_temp, patch = prod_temp)
plotAtlas(prod_atlas)
close3d()

###---------------------------------------------------------------------------

## STEP 6.5 Place the patch onto all specimens.
# Deform the template over each mesh surface. The inflate parameter can be adjusted if there are placement errors observed in the checking stage
prod_patched <- placePatch(atlas = prod_atlas, dat.array = prod_outline, path = "analysis/data/raw_data/product_data", inflate = 5, fileext = ".ply")

###---------------------------------------------------------------------------

## STEP 6.6 Check patch slms have deformed correctly on each specimen (this is time consuming but essential to check the deformation process). If there are placement errors, try increasing the inflate parameter, or repeat previous steps excluding problem specimens.
checkLM(prod_patched, atlas=prod_atlas)
close3d()

# Save the patched slm data
save(prod_patched, file = "analysis/data/derived_data/prod_patched.RData")

## Dataset is ready for GM comparison of end-products. The next step is required to extract the product ventral for comparison with the final preferential scar on cores.

###---------------------------------------------------------------------------

## STEP 6.7 Extract outline of product ventral
# Create a new empty data array. The percussion point (1) and dorsal platform points (49-58) are not needed, reducing the number of lms to 52.
ventr_lm <- array(NA, dim = c(52, 3, 178))

# Reconfigure original product landmarks in the new array to only include the ventral outline and reordering them
ventr_lm[1:40, , ] <- prod_landmarks[4:43, , ] # (outline)
ventr_lm[41, , ] <- prod_landmarks[2, , ]   # (R platform end)
ventr_lm[42:46, , ] <- prod_landmarks[59:63, , ]   # (platform ventr)
ventr_lm[47:51, , ] <- prod_landmarks[44:48, , ]  # (platform ventr)
ventr_lm[52, , ] <- prod_landmarks[3, , ] # (L platform end)

# Create a new data array for resampled ventral outline semilandmarks (37 slms to match the core pref scar)
prod_ventr <- array(0, dim=c(37,3,178))

# Redistribute and resample equally-spaced semilandmarks on all specimens along two curves (5 platform and 30 edge points)
for(i in 1:dim(ventr_lm)[3]) {
  startA <- ventr_lm[52,,i]
  startB <- ventr_lm[41,,i]
  curveA <- digit.curves(startA,rbind(ventr_lm[1:40,,i], startB), nPoints = 30, closed = F)
  curveB <- digit.curves(startB, rbind(ventr_lm[42:51,,i], startA), nPoints = 5, closed = F)
  comb <- rbind(curveA[2:nrow(curveA),], curveB[2:nrow(curveB),])
  prod_ventr[,,i] <- comb
}

# View equally spaced slms on mesh model
shade3d(product1, color = "gray")
spheres3d(prod_ventr[31,,3],col=1,radius=1)
spheres3d(prod_ventr[37,,3],col=1,radius=1)
spheres3d(prod_ventr[1:30,,3],col=3,radius=0.7)
spheres3d(prod_ventr[32:36,,3],col=4,radius=0.7)
#rgl.snapshot('analysis/figures/prodventr_landmarks_example.png', fmt = 'png')
close3d()

# Define fixed LMs (that do not slide) and semilandmarks on outline curves to relax landmarks
fixPv <- c(31,37)
outlineA <- c(1:30)
outlineB <- c(32:36)

outlinesPv <- list(outlineA, outlineB)

dimnames(prod_ventr)[[3]] <- names

# Slide semilandmarks along curves (this step takes some time and processing power for the large sample size)
slide_curvesPv <- slider3d(prod_ventr, SMvector = fixPv, deselect = TRUE, outlines = outlinesPv,
                          sur.path = "analysis/data/raw_data/product_data", sur.type = "ply", iterations = 3)

# View slid slms on mesh model
shade3d(product1, color = "gray")
spheres3d(slide_curvesPv$dataslide[fixPv,,3],col=4,radius=1)
spheres3d(slide_curvesPv$dataslide[,,3],col=3,radius=0.7)
close3d()

# Rename and save
prod_voutline <- slide_curvesPv$dataslide
save(prod_voutline, file = "analysis/data/derived_data/prod_voutline.RData")

###---------------------------------------------------------------------------

## STEP 6.8 Append attribute data to specimens

# Import attribute data
pts <- read.csv("File_list_prods.csv")
pts <- dplyr::arrange(pts, File_name)
Sample <- factor(pts$Assemblage)
Platf_type <- factor(pts$Platf_type)
Dist_term <- factor(pts$Dist_term)
Flaking_axis <- factor(pts$Flaking_axis)
Flake_scars <- factor(pts$Flake_scars)

###---------------------------------------------------------------------------

## Product data are now ready for GM analyses

## R data files can be loaded as:
load("analysis/data/derived_data/prod_patched.RData") # product 3D surfaces
load("analysis/data/derived_data/prod_outline.RData") # product ventral outline only (for comparison with core pref scar)

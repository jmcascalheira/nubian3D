
## STEP 7. R SCRIPT FOR ALLOMETRY ANALYSES

# Ensure all required libraries are loaded
library(geomorph)
library(ggplot2)
library(dplyr)

###---------------------------------------------------------------------------

# Load R data files of core landmark configurations generated in Step 5:
load("analysis/data/derived_data/allpatched.RData")
load("analysis/data/derived_data/alloutline.RData")
load("analysis/data/derived_data/allprefs.RData")

# Import the attribute data for cores and read attributes as Factors
cores <- read.csv("analysis/data/derived_data/File_list_cores.csv", header=TRUE, stringsAsFactors = TRUE)

###---------------------------------------------------------------------------

## -- CORE SURFACES --
## STEP 7.1 For all core surfaces ("allpatched"), perform GPA and create new geomorph dataframe
# Perform Generalised Procrustes Analysis (GPA)
patched_gpa<-gpagen(allpatched, print.progress = FALSE)

# Create a new geomorph dataframe (.gdf) from the GPA output
allpatch.gdf <- geomorph.data.frame(patched_gpa)

# Append the attribute/factor data to the gdf
allpatch.gdf <- append(allpatch.gdf, cores)

# Creates new column in the gdf with natural log of the centroid size
allpatch.gdf$lnCS <- log(allpatch.gdf$Csize)

###---------------------------------------------------------------------------

## STEP 7.2 Compare centroid sizes
# Create new table containing core attribute data and centroid size from the GPA data
coresize <- mutate(cores, Csize = patched_gpa$Csize)

# Summarise descriptive statistics (mean, min, max, sd, Coefficient of Variation (CV))
Csize_stats_Assemblage <- coresize %>%
  group_by(Assemblage) %>%
  summarise(
    mean_Csize = mean(Csize, na.rm = TRUE),
    min_Csize = min(Csize, na.rm = TRUE),
    max_Csize = max(Csize, na.rm = TRUE),
    sd_Csize = sd(Csize, na.rm = TRUE),
    CV_Csize = (sd_Csize / mean_Csize) * 100
  )
print(Csize_stats_Assemblage)

Csize_stats_Area <- coresize %>%
  group_by(Area) %>%
  summarise(
    mean_Csize = mean(Csize, na.rm = TRUE),
    min_Csize = min(Csize, na.rm = TRUE),
    max_Csize = max(Csize, na.rm = TRUE),
    sd_Csize = sd(Csize, na.rm = TRUE),
    CV_Csize = (sd_Csize / mean_Csize) * 100
  )
print(Csize_stats_Area)


# T-test to compare centroid size for both regions
t_testArea <- t.test(Csize ~ Area, data = coresize)
print(t_testArea)


# Plot centroid size boxplot by assemblage
ggplot(coresize, aes(x = Assemblage, y = Csize, fill = Area)) +
  geom_jitter(position = position_jitter (width = 0.4), size = 2, alpha = 0.6, aes(color = Area), show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE, alpha = 0.7) +
  labs(x = "Region", y = "Centroid size") +
  scale_fill_manual(values = c(NK = "royalblue", TH = "indianred1")) +
  scale_color_manual(values = c(NK = "royalblue", TH = "indianred1")) +
  theme_minimal()


# Create new table containing pref scar attribute data and centroid size from the GPA data
pref_gpa<-gpagen(allprefs, print.progress = FALSE)
prefsize <- mutate(cores, Csize = pref_gpa$Csize)

# Summarise descriptive statistics (mean, min, max)
Csize_pref_Assemblage <- prefsize %>%
  group_by(Assemblage) %>%
  summarise(
    mean_Csize = mean(Csize, na.rm = TRUE),
    min_Csize = min(Csize, na.rm = TRUE),
    max_Csize = max(Csize, na.rm = TRUE),
    sd_Csize = sd(Csize, na.rm = TRUE),
    CV_Csize = (sd_Csize / mean_Csize) * 100
  )
print(Csize_pref_Assemblage)

Csize_pref_Area <- prefsize %>%
  group_by(Area) %>%
  summarise(
    mean_Csize = mean(Csize, na.rm = TRUE),
    min_Csize = min(Csize, na.rm = TRUE),
    max_Csize = max(Csize, na.rm = TRUE),
    sd_Csize = sd(Csize, na.rm = TRUE),
    CV_Csize = (sd_Csize / mean_Csize) * 100
  )
print(Csize_pref_Area)

###---------------------------------------------------------------------------

## STEP 7.3 Perform the allometry analyses for core surfaces
# Perform Procrustes ANOVA to regress co-ordinates against log centroid size
patchCS <- procD.lm(coords~ log(Csize), data=allpatch.gdf, iter=1000, RRPP=TRUE, print.progress = FALSE)

# View ANOVA results
summary(patchCS)

# Define colors and shapes for NK and TH
color <- c(NK = "royalblue", TH = "indianred1")
shape <- c(NK = 16, TH = 17)

# Assign colors and shapes based on Area in the .gdf
allpatch.gdf$Area <- as.factor(allpatch.gdf$Area)
color <- color[allpatch.gdf$Area]
shape <- shape[allpatch.gdf$Area]

# Plot regression score against centroid size (i.e. show allometry)
RegPatch <- plotAllometry(patchCS, size=allpatch.gdf$Csize,  method = "RegScore", col=color, pch=shape, bg=allpatch.gdf$Area, xlab="Ln Core Surface Centroid Size")

# Select specific data points to view deformation grid of corresponding specimens in rgl viewer. These can be output as .png
picknplot.shape(RegPatch)

###---------------------------------------------------------------------------

## -- CORE OUTLINES --
## STEP 7.4 Repeat the GPA and Procrustes ANOVA process for all core outlines ("alloutline")
# Perform GPA and create new geomorph dataframe with attribute and log centroid data
outline_gpa <- gpagen(alloutline, print.progress = FALSE)
alloutline.gdf <- geomorph.data.frame(outline_gpa)
alloutline.gdf <- append(alloutline.gdf, cores)
alloutline.gdf$lnCS <- log(alloutline.gdf$Csize)

# Perform Procrustes ANOVA and generate plots for all core outlines
outlineCS <- procD.lm(coords~ log(Csize), data=alloutline.gdf, iter=1000, RRPP=TRUE, print.progress = FALSE)

# View ANOVA results
summary(outlineCS)

# Assign colors and shapes based on Area in the .gdf
alloutline.gdf$Area <- as.factor(alloutline.gdf$Area)
color <- colsArea[alloutline.gdf$Area]
shape <- shapesArea[alloutline.gdf$Area]

# Plot Allometry
RegOutline <- plotAllometry(outlineCS, size=alloutline.gdf$Csize,  method = "RegScore", col=color, pch=shape, bg=alloutline.gdf$Area, xlab="Ln Core Outline Centroid Size")

# View specific specimens
picknplot.shape(RegOutline)

###---------------------------------------------------------------------------

## -- PREF SCARS --
## STEP 7.5 Repeat the GPA process for all core preferential scars ("allprefs")
# Perform GPA and create new geomorph dataframe with attribute and log centroid data
pref_gpa <-gpagen(allprefs, print.progress = FALSE)
allprefs.gdf <- geomorph.data.frame(pref_gpa)
allprefs.gdf <- append(allprefs.gdf, cores)
allprefs.gdf$lnCS <- log(allprefs.gdf$Csize)

# Perform Procrustes ANOVA and generate plots for all core preferential scars
prefCS <- procD.lm(coords~ log(Csize), data=allprefs.gdf, iter=1000, RRPP=TRUE, print.progress = FALSE)

# View ANOVA results
summary(prefCS)

# Set plot aesthetics
allprefs.gdf$Area <- as.factor(allprefs.gdf$Area)
color <- colsArea[allprefs.gdf$Area]
shape <- shapesArea[allprefs.gdf$Area]

# Plot allometry
RegPref <- plotAllometry(prefCS, size=allprefs.gdf$Csize,  method = "RegScore", col=color, pch=shape, bg=allprefs.gdf$Area, xlab="Ln Pref Scar Centroid Size")

# View specific specimens
picknplot.shape(RegPref)

###---------------------------------------------------------------------------

## STEP 7.6 (Optional) Plot allometry plots in ggplot2 to show regression line. Note that geomorph v4 has a function to convert plots to ggplot objects where compatible
# Core surfaces ("RegPatch")
datapatch <- data.frame(lncsize = log(allpatch.gdf$Csize), regsc = RegPatch$RegScore)

ggplot(datapatch, aes(x = lncsize, y = regsc)) +
  geom_point(aes(col = allpatch.gdf$Area, shape = allpatch.gdf$Area)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE) +
  labs(x = "Ln Core Surface Centroid Size", y = "Regression Score", color = "Region", shape = "Region") +
  scale_color_manual(values = colsArea) +
  scale_shape_manual(values = shapesArea)+
  geom_hline(yintercept = -0.25, linetype = "solid", color = "darkgrey") + # Horizontal axis line
  geom_vline(xintercept = 5.5, linetype = "solid", color = "darkgrey") + # Vertical axis line
  theme_minimal()+
  ggtitle("Core surfaces: allometry plot")

# Calculate slope of regression lines on plot
regline <- lm(regsc ~ lncsize, data = datapatch)
slope <- coefficients(regline)[2]
print(slope)


# Core outlines ("RegOutline")
dataoutl <- data.frame(lncsize = log(alloutline.gdf$Csize), regsc = RegOutline$RegScore)

ggplot(dataoutl, aes(x = lncsize, y = regsc)) +
  geom_point(aes(col = allpatch.gdf$Area, shape = allpatch.gdf$Area)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE) +
  labs(x = "Ln Core Outline Centroid Size", y = "Regression Score", color = "Region", shape = "Region") +
  scale_color_manual(values = colsArea) +
  scale_shape_manual(values = shapesArea)+
  geom_hline(yintercept = -0.2, linetype = "solid", color = "darkgrey") + # Horizontal axis line
  geom_vline(xintercept = 4.5, linetype = "solid", color = "darkgrey") + # Vertical axis line
  theme_minimal()+
  ggtitle("Core outlines: allometry plot")

# Calculate slope of regression lines on plot
regline <- lm(regsc ~ lncsize, data = dataoutl)
slope <- coefficients(regline)[2]
print(slope)


# Core pref scar ("RegPref")
datapref <- data.frame(lncsize = log(allprefs.gdf$Csize), regsc = RegPref$RegScore)

ggplot(datapref, aes(x = lncsize, y = regsc)) +
  geom_point(aes(col = allpatch.gdf$Area, shape = allpatch.gdf$Area)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5, fullrange = TRUE) +
  labs(x = "Ln Pref Scar Centroid Size", y = "Regression Score", color = "Region", shape = "Region") +
  scale_color_manual(values = colsArea) +
  scale_shape_manual(values = shapesArea)+
  geom_hline(yintercept = -0.4, linetype = "solid", color = "darkgrey") + # Horizontal axis line
  geom_vline(xintercept = 4, linetype = "solid", color = "darkgrey") + # Vertical axis line
  theme_minimal()+
  ggtitle("Pref scars: allometry plot")

# Calculate slope of regression lines on plot
regline <- lm(regsc ~ lncsize, data = datapref)
slope <- coefficients(regline)[2]
print(slope)

###---------------------------------------------------------------------------

## STEP 7.7  Further investigate regional differences by conducting GPA and Procrustes ANOVA for each region separately
## -- CORE SURFACES --
# Load NK and TH data
load("patchedNK.RData")
load("patchedTH.RData")

# Perform Generalised Procrustes Analysis (GPA)
patchedNK_gpa<-gpagen(patchedNK, print.progress = FALSE)
patchedTH_gpa<-gpagen(patchedTH, print.progress = FALSE)

# Create a new geomorph dataframe (.gdf) from the GPA output
NKpatch.gdf <- geomorph.data.frame(patchedNK_gpa)
THpatch.gdf <- geomorph.data.frame(patchedTH_gpa)

# Perform Generalised Procrustes Analysis (GPA)
prefNK_gpa<-gpagen(prefNK, print.progress = FALSE)
prefTH_gpa<-gpagen(prefTH, print.progress = FALSE)

# Perform Procrustes ANOVA for NK core surfaces
NK_CS <- procD.lm(coords ~ log(Csize), data = NKpatch.gdf, iter = 1000, RRPP = TRUE, print.progress = FALSE)
summary(NK_CS)

# Perform Procrustes ANOVA for TH core surfaces
TH_CS <- procD.lm(coords ~ log(Csize), data = THpatch.gdf, iter = 1000, RRPP = TRUE, print.progress = FALSE)
summary(TH_CS)


## -- PREF SCARS --
# Load NK and TH data
load("prefNK.RData")
load("prefTH.RData")

# Perform Generalised Procrustes Analysis (GPA)
prefNK_gpa<-gpagen(prefNK, print.progress = FALSE)
prefTH_gpa<-gpagen(prefTH, print.progress = FALSE)

# Create a new geomorph dataframe (.gdf) from the GPA output
NKpref.gdf <- geomorph.data.frame(prefNK_gpa)
THpref.gdf <- geomorph.data.frame(prefTH_gpa)

# Perform Procrustes ANOVA for NK preferential scars
NKpr_CS <- procD.lm(coords ~ log(Csize), data = NKpref.gdf, iter = 1000, RRPP = TRUE, print.progress = FALSE)
summary(NKpr_CS)

# Perform Procrustes ANOVA for TH preferential scars
THpr_CS <- procD.lm(coords ~ log(Csize), data = THpref.gdf, iter = 1000, RRPP = TRUE, print.progress = FALSE)
summary(THpr_CS)

## PC PLOTS SUPPLEMENTARY TO STEP 8
# These plots are presented in figures/supplementary figures but not included in the main Step 8 script

# Set colour/shape aesthetics for the NK and TH samples for the plots
colsArea <- c("NK" = "royalblue", "TH" = "indianred1")
shapesArea <- c("NK" = 16, "TH" = 17)

# Set colour aesthetics for preparation (scar direction) and Nubian core Type
colsScars <- c("Distal" = "#FC8D62", "Bilateral" = "#66C2A5", "Distal-bilateral" = "#A6D854", "Distal-lateral-L" = "#E78AC3", "Distal-lateral-R" = "#8DA0CB")
colsProdscar <- c("Bidirectional" = "#FC8D62", "Bilateral" = "#66C2A5", "Unidirectional" = "#FFD92F", "Unidirectional convergent" = "#FFD92F", "Centripetal" = "#A6D854", "Crossed L" = "#E78AC3", "Crossed R" = "#8DA0CB", "Indeterminate" = "grey")
colsType <- c("T1" = "#FC8D62", "T2" = "#66C2A5", "T1/2" = "#8DA0CB")

# Set colour/shape aesthetics for assemblage
colsAssembl <- c("NK1_M" = "#80B1D3", "NK1_U" = "#BEBADA", "NK3" = "#8DD3C7", "TH584" = "#FB8072", "TH571" = "#FDB462")
shapeAssembl <- c("NK1_M" = 16, "NK1_U" = 17, "NK3" = 18)

###---------------------------------------------------------------------------

## -- CORE SURFACES -- STEP 8.4
# Aesthetics set for Assemblage (Figure S1)
#PC1 and 2
coord_allpatch %>%
  ggplot(aes(x = Comp1, y = Comp2), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Area)) +
  scale_color_manual(values = colsAssembl) +
  scale_fill_manual(values = colsAssembl) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Assemblage") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")

#PC3 and 4
coord_allpatch %>%
  ggplot(aes(x = Comp3, y = Comp4), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Area)) +
  scale_color_manual(values = colsAssembl) +
  scale_fill_manual(values = colsAssembl) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Assemblage") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC3", y = "PC4") +
  theme(legend.position = "right")

#PC5 and 6
coord_allpatch %>%
  ggplot(aes(x = Comp5, y = Comp6), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Area)) +
  scale_color_manual(values = colsAssembl) +
  scale_fill_manual(values = colsAssembl) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Assemblage") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC5", y = "PC6") +
  theme(legend.position = "right")


# Aesthetics set for Core Type (Figure S2)
#PC1 and 2
coord_allpatch %>%
  ggplot(aes(x = Comp1, y = Comp2), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Type, shape = Area)) +
  scale_color_manual(values = colsType) +
  scale_fill_manual(values = colsType) +
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Type") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")

#PC3 and 4
coord_allpatch %>%
  ggplot(aes(x = Comp3, y = Comp4), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Type, shape = Area)) +
  scale_color_manual(values = colsType) +
  scale_fill_manual(values = colsType) +
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Type") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC3", y = "PC4") +
  theme(legend.position = "right")


# Aesthetics set for preparation/scar pattern (plots not presented in paper)
#PC1 and 2
coord_allpatch %>%
  ggplot(aes(x = Comp1, y = Comp2), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Preparation, shape = Area)) +
  scale_color_manual(values = colsScars) +
  scale_fill_manual(values = colsScars) +
  stat_ellipse(geom = "polygon", aes(fill = Preparation), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Preparation") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")

#PC3 and 4
coord_allpatch %>%
  ggplot(aes(x = Comp3, y = Comp4), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Preparation, shape = Area)) +
  scale_color_manual(values = colsScars) +
  scale_fill_manual(values = colsScars) +
  stat_ellipse(geom = "polygon", aes(fill = Preparation), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Preparation") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC3", y = "PC4") +
  theme(legend.position = "right")

#PC5 and 6
coord_allpatch %>%
  ggplot(aes(x = Comp5, y = Comp6), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Preparation, shape = Area)) +
  scale_color_manual(values = colsScars) +
  scale_fill_manual(values = colsScars) +
  stat_ellipse(geom = "polygon", aes(fill = Preparation), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Preparation") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC5", y = "PC6") +
  theme(legend.position = "right")

###---------------------------------------------------------------------------

## -- PREF SCARS -- STEP 8.6
# Aesthetics set for area/region (Figure 7)
#PC1 and 3
coord_allprefs %>%
  ggplot(aes(x = Comp1, y = Comp3, shape = Area)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Area)) +
  scale_color_manual(values = colsArea) +
  scale_fill_manual(values = colsArea) +
  stat_ellipse(geom = "polygon", aes(fill = Area), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Area") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

#PC2 and 4 (asymmetry)
coord_allprefs %>%
  ggplot(aes(x = Comp2, y = Comp4, shape = Area)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Area)) +
  scale_color_manual(values = colsArea) +
  scale_fill_manual(values = colsArea) +
  stat_ellipse(geom = "polygon", aes(fill = Area), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Area") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC2", y = "PC4") +
  theme(legend.position = "right")


# Aesthetics set for Type (plots not presented in paper)
#PC1 and 3
coord_allprefs %>%
  ggplot(aes(x = Comp1, y = Comp2), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Type, shape = Area)) +
  scale_color_manual(values = colsType) +
  scale_fill_manual(values = colsType) +
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Type") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

#PC2 and 4
coord_allprefs %>%
  ggplot(aes(x = Comp3, y = Comp4), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Type, shape = Area)) +
  scale_color_manual(values = colsType) +
  scale_fill_manual(values = colsType) +
  stat_ellipse(geom = "polygon", aes(fill = Type), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Type") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC2", y = "PC4") +
  theme(legend.position = "right")


# Aesthetics set for preparation (plots not presented in paper)
#PC1 and 3
coord_allprefs %>%
  ggplot(aes(x = Comp1, y = Comp3), shape = Area) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Preparation, shape = Area)) +
  scale_color_manual(values = colsScars) +
  scale_fill_manual(values = colsScars) +
  stat_ellipse(geom = "polygon", aes(fill = Preparation), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Preparation") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

###---------------------------------------------------------------------------

## -- NK PRODUCTS -- STEP 8.7
# Aesthetics set for Assemblage (Figure 8)
#PC 1 and 3
coord_prods %>%
  ggplot(aes(x = Comp1, y = Comp3), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Assemblage)) +
  #scale_color_manual(values = colsAssembl) +
  #scale_fill_manual(values = colsAssembl) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Assemblage") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

#PC 2 and 4
coord_prods %>%
  ggplot(aes(x = Comp2, y = Comp4), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Assemblage)) +
  #scale_color_manual(values = colsAssembl) +
  #scale_fill_manual(values = colsAssembl) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  labs(fill = "Assemblage") +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC2", y = "PC4") +
  theme(legend.position = "right")


# Aesthetics set for scar pattern (plots not presented in paper)
#PC 1 and 3
coord_prods %>%
  ggplot(aes(x = Comp1, y = Comp3), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Flake_scars, shape = Assemblage)) +
  scale_color_manual(values = colsProdscar) +
  scale_fill_manual(values = colsProdscar) +
  stat_ellipse(geom = "polygon", aes(fill = Flake_scars), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

#PC 2 and 4
coord_prods %>%
  ggplot(aes(x = Comp2, y = Comp4), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Flake_scars, shape = Assemblage)) +
  scale_color_manual(values = colsProdscar) +
  scale_fill_manual(values = colsProdscar) +
  stat_ellipse(geom = "polygon", aes(fill = Flake_scars), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC2", y = "PC4") +
  theme(legend.position = "right")

#PC 5 and 6
coord_prods %>%
  ggplot(aes(x = Comp5, y = Comp6), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Flake_scars, shape = Assemblage)) +
  scale_color_manual(values = colsProdscar) +
  scale_fill_manual(values = colsProdscar) +
  stat_ellipse(geom = "polygon", aes(fill = Flake_scars), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC5", y = "PC6") +
  theme(legend.position = "right")


# -- PREF SCAR AND PRODUCT VENTRAL --
#PC 1 and 3 by artefact (Figure S3)
coord_prod_scar %>%
  ggplot(aes(x = Comp1, y = Comp3), shape = Assemblage) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Artefact, shape = Assemblage)) +
  stat_ellipse(geom = "polygon", aes(fill = Artefact), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

#PC 1 and 3 by assemblage (Figure S3)
coord_prod_scar %>%
  ggplot(aes(x = Comp1, y = Comp3), shape = Artefact) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = Assemblage, shape = Artefact)) +
  stat_ellipse(geom = "polygon", aes(fill = Assemblage), alpha = 0.2, show.legend = TRUE, level = 0.95) +
  theme_minimal() +
  ggtitle("") +
  labs(x = "PC1", y = "PC3") +
  theme(legend.position = "right")

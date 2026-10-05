############# SPATIAL
import_summary <- tibble(
  "Method" = c("cellchat", "lianaP_morans", "mistyR", "NICHES", "spatialDM", "stlearn", "seurat_wilcoxon", "cellphonedbv5"),
  "How it computes interaction" = c(
    "Law of mass action, multiplying the average expression values of all ligand/receptor subunits, cofactors, and complexes. For spatial context, the expression product is explicitly modulated by a continuous distance decay function",
    "It employs global spatial bivariate autocorrelation statistic (Moran's I) across entire tissue. It also allows for local neighborhood-weighted cosine similarity to be calculated (not used here)",
    "For each marker of interest, mistyR builds prediction models (linear model used in this benchmark but in general is algorithm-agnostic framework) that relate its expression to other markers in different spatial contexts",
    "NICHES preserves single-cell resolution and computes spatial interaction for each pairs of single-cells within radius. For each edge, interaction score is calculated by cross-product of the normalized ligand and receptor expression in corresponding cells",
    "SpatialDM identifies significant interactions by detecting spatially co-expressed ligand-receptor pairs through bivariate Moran's I (Moran's R)",
    "Interaction scores are calculated for each spot based on expression of ligand and the expression of receptor in all surrounding spots",
    "Non spatially aware method",
    "Computes interaction based on average expression across celltypes defined in the microenvironment file provided by the user"
  ),
  "Statistics" = c(
    "Permutation of cell labels",
    "Permutation of gene expression vectors",
    "mistyR retrieves feature importance as opposed to pvalus. The importances depend on prediction model used. For linear model feature importances are defines as t-values of the model coefficients calculated via Ordinary Least Squares (OLS)",
    "Significance of interaction is determined using FindAllMarkers function from Seurat",
    "Permutation by shuffling cell IDs",
    "stlearn uses 2-level permutation test. 
    LR significance testing: generates a random background signal of non-interaction genes. P-values for each spot and LR pair are derived based on proportion of the background scores that had a score greater than the LRscore
    Celltype specific testing: For specific satial locations obtained from previous testing, stlearn permutes celltype information to assess if specific celltypes are over-represented",
    "Wilcoxon test",
    "Permutation of cell labels"
  )
)

# Render block fixed for Markdown rendering compatibility
import_summary_table <- import_summary %>%
  gt() %>%  # FIXED: Removed groupname_col = "Method" so Markdown stays active
  fmt_markdown(columns = `How it computes interaction`) %>% # REQUIRED: Tells gt to process the <br> tag
  tab_header(
    title = md("**Table 2. Review of spatial computational methods for cell-cell communication**")
  ) %>%
  # Journal Style Sheet Rules (Nature Style Elements)
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1.5),
    table_body.border.bottom.color = "black",
    table_body.border.bottom.width = px(1.5)
  ) %>%
  cols_align(
    align = "left", # Better readability for long character descriptions
    columns = `How it computes interaction`
  ) %>%
  cols_align(
    align = "center",
    columns = c(Method, Statistics)
  )

# Display the table in your Quarto preview
import_summary_table


import_summary_table %>%
  gtsave(
    filename = "Table2_spatial_Cell_Communication_Summary.pdf"
  )

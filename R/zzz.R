.onAttach <- function(libname, pkgname) {
  
  # Get package version dynamically
  version <- utils::packageVersion(pkgname)
  
  # 1. Brief Description
  packageStartupMessage(
    "Loading spatialCooccur v", version, ": An R package for analyzing spatial co-occurrence."
  )
  
  # 2. Citation Information
  citation_text <- "Inamo J, et al. (2026). Spatial transcriptomics reveals immune-stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis. JCI Insight. 11(1):e198074. doi:10.1172/jci.insight.198074"
  packageStartupMessage(
    "To cite this package in publications, please use:\n",
    strwrap(citation_text, prefix = "  ") # strwrap formats long text nicely
  )

  # 3. Developer and Contact Information
  packageStartupMessage(
    "Developed by: Jun Inamo <juninamo@keio.jp>"
  )
}
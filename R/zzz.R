.onAttach <- function(libname, pkgname) {
  
  # Get package version dynamically
  version <- utils::packageVersion(pkgname)
  
  # 1. Brief Description
  packageStartupMessage(
    "Loading spatialCooccur v", version, ": An R package for analyzing spatial co-occurrence."
  )
  
  # 2. Citation Information
  citation_text <- "Inamo J, et al. (2025). Subcellular spatial transcriptomics reveals immune-stromal crosstalk within the synovium of patients with juvenile idiopathic arthritis. medRxiv. doi:10.1101/2025.08.05.25332835"
  packageStartupMessage(
    "To cite this package in publications, please use:\n",
    strwrap(citation_text, prefix = "  ") # strwrap formats long text nicely
  )
  
  # 3. Developer and Contact Information
  packageStartupMessage(
    "Developed by: Jun Inamo <jun.inamo@cuanschutz.edu>"
  )
}
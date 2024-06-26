% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convPairs.R
\name{convPairs}
\alias{convPairs}
\title{Converts colPairs object into an abundance matrix based on number of nearby
interactions for every cell type.}
\usage{
convPairs(cells, colPair, cellType = "cellType", imageID = "imageID")
}
\arguments{
\item{cells}{A SingleCellExperiment that contains objects in the colPairs slot.}

\item{colPair}{The name of the object in the colPairs slot for which the dataframe is
constructed from.}

\item{cellType}{The cell type if using SingleCellExperiment.}

\item{imageID}{The image ID if using SingleCellExperiment.}
}
\value{
Matrix of abundances
}
\description{
Converts colPairs object into an abundance matrix based on number of nearby
interactions for every cell type.
}
\examples{
data("diabetesData")
images <- c("A09", "A11", "A16", "A17")
diabetesData <- diabetesData[
  , SummarizedExperiment::colData(diabetesData)$imageID \%in\% images
]

diabetesData_SPE <- SpatialExperiment::SpatialExperiment(diabetesData,
  colData = SummarizedExperiment::colData(diabetesData)
)
SpatialExperiment::spatialCoords(diabetesData_SPE) <- data.frame(
  SummarizedExperiment::colData(diabetesData_SPE)$x,
  SummarizedExperiment::colData(diabetesData_SPE)$y
) |>
  as.matrix()

SpatialExperiment::spatialCoordsNames(diabetesData_SPE) <- c("x", "y")

diabetesData_SPE <- imcRtools::buildSpatialGraph(diabetesData_SPE,
  img_id = "imageID",
  type = "knn",
  k = 20,
  coords = c("x", "y")
)

pairAbundances <- convPairs(diabetesData_SPE,
  colPair = "knn_interaction_graph"
)

}

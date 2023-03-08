#' Metadata information
#'
#' This data frame contains information about the samples used in the study.
#'
#' @format A data frame with 14 variables:
#' \describe{
#' \item{BW}{Body weight of the sample in grams.}
#' \item{Length}{Length of the sample in centimeters.}
#' \item{Sex}{Sex of the sample.}
#' \item{Gut_struc}{Structural condition of the gut of the sample.}
#' \item{Gut_light}{Light coloration of the gut of the sample.}
#' \item{Anus_bleed}{Presence or absence of anus bleeding in the sample.}
#' \item{Gut_bleed}{Presence or absence of gut bleeding in the sample.}
#' \item{Fat_Score}{Score of fat deposition in the sample.}
#' \item{Heart_Sco}{Score of heart condition in the sample.}
#' \item{softliver}{Presence or absence of soft liver in the sample.}
#' \item{patchy_liver}{Presence or absence of patchy liver in the sample.}
#' \item{Liver_Sco}{Score of liver condition in the sample.}
#' \item{HSI}{Hepatosomatic index of the sample.}
#' \item{CSI}{Cardiosomatic index of the sample.}
#' }
#'
#' @source Information was collected from the study samples.
#'
#' @examples
#' # View the first few rows of the metadata
#' head(metadata)
#'
#' @export
metadata <- data.frame(
  BW = numeric(),
  Length = numeric(),
  Sex = character(),
  Gut_struc = character(),
  Gut_light = character(),
  Anus_bleed = character(),
  Gut_bleed = character(),
  Fat_Score = numeric(),
  Heart_Sco = numeric(),
  softliver = character(),
  patchy_liver = character(),
  Liver_Sco = numeric(),
  HSI = numeric(),
  CSI = numeric()
)

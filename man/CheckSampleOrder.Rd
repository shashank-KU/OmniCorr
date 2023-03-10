% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CheckSampleOrder.R
\name{CheckSampleOrder}
\alias{CheckSampleOrder}
\title{Check if two data frames have matching sample names and order}
\usage{
CheckSampleOrder(df1, df2)
}
\arguments{
\item{df1}{A data frame with data, where rows represent samples}

\item{df2}{A data frame with data, where rows represent samples}
}
\value{
A list with both df1 and df2 data in the same order and with matching rownames
}
\description{
This function checks if two data frames have matching rownames and row order. If the rownames
match but are in a different order, the function will reorder both data frames to match.
}
\examples{
data1 <- data.frame(A=c(1,2,3), B=c(4,5,6), row.names=c("S1", "S2", "S3"))
data2 <- data.frame(X=c(7,8,9), Y=c(10,11,12), row.names=c("S2", "S1", "S3"))
df_list <- CheckSampleOrder(data1, data2)
data1_reordered <- df_list[[1]]
data2_reordered <- df_list[[2]]
}

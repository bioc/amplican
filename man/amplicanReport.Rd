% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amplicanReport.R
\name{amplicanReport}
\alias{amplicanReport}
\title{Prepare reports as .Rmd files.}
\usage{
amplicanReport(results_folder, levels = c("id", "barcode", "group", "guide",
  "amplicon", "summary"), report_files = c("report_id", "report_barcode",
  "report_group", "report_guide", "report_amplicon", "report_summary"))
}
\arguments{
\item{results_folder}{(string) Folder containing results from the
\code{\link{amplicanAlign}} function, do not change names of the files.}

\item{levels}{(vector) Possible values are: 'id', 'barcode', 'group',
'guide', 'amplicon', 'summary'. You can also input more than one value
eg. c('id', 'barcode') will create two separate reports for each level.}

\item{report_files}{(vector) You can supply your own names of the files.
For each of the levels there has to be one file name. Files are created
in current working directory by default.}
}
\value{
NULL All results are written into current dir with specified names.
}
\description{
amplicanReport takes a configuration file, fastq reads and output directory
to prepare summary as editable .Rmd file. You can specify whether you want
to make summaries based on ID, Barcode, Group or even guideRNA and Amplicon.
}
\examples{
results_folder <- system.file("extdata", "results", package = "amplican") #output folder
amplicanReport(results_folder, report_files = paste0(paste0(results_folder, "/reports"),
                                                     c("/report_id",
                                                       "/report_barcode",
                                                       "/report_group",
                                                       "/report_guide",
                                                       "/report_amplicon",
                                                       "/report_summary")))

}
\seealso{
Other analysis steps: \code{\link{amplicanAlign}},
  \code{\link{amplicanPipeline}}
}

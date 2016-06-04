\name{IFNB_Data}
\alias{IFNB_Data}
\docType{data}
\title{IFNB Gene Time Course Dataset}
\description{
Gene expression time-course data from Baranzini et al (2005), for patients with multiple sclerosis after treatment with IFNBeta.  Patients are classified as good or bad responders. Gene expression was measuredfor 76 genes of interest before treatment (baseline) and at 6 follow-up time points over the nexttwo years (3 months, 6 months, 9 months, 12 months, 18 months, 24 months).  The original data has been log-normalized and mean-imputed.          
}
\usage{data(IFNB_Data)}
\format{
\itemize{
\item DataArray: Time course array: 53 patients X 7 time points X 76 genes.  
\item Class: Vector giving patient classification as good or poor responders.  
}}
   
\references{
Baranzini, S.E, Mousavi, P., Rio, J., et al. (2005). Transcription-based prediction of responseto IFNBeta using supervised computational methods. PLoS Biology 3(1), e2.}

\keyword{datasets}

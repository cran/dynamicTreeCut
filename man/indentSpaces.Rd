\name{indentSpaces}
\alias{indentSpaces}
\title{Spaces for indented output. }
\description{Returns a character string containing two times \code{indent} spaces.
}
\usage{
indentSpaces(indent = 0)
}
\arguments{
  \item{indent}{Desired level of indentation. The number of returned spaces will be twice this
argument.}
}
\value{A character string containing spaces, of length twice \code{indent}. }
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\examples{
spaces = indentSpaces(0);
print(paste(spaces, "This output is not indented..."));
spaces = indentSpaces(1);
print(paste(spaces, "...while this one is."))
}

\keyword{print}

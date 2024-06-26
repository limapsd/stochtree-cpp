% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{CppRNG}
\alias{CppRNG}
\title{Class that wraps a C++ random number generator (for reproducibility)}
\description{
Persists a C++ random number generator throughout an R session to
ensure reproducibility from a given random seed. If no seed is provided,
the C++ random number generator is initialized using \code{std::random_device}.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{rng_ptr}}{External pointer to a C++ std::mt19937 class}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-CppRNG-new}{\code{CppRNG$new()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-CppRNG-new"></a>}}
\if{latex}{\out{\hypertarget{method-CppRNG-new}{}}}
\subsection{Method \code{new()}}{
Create a new CppRNG object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CppRNG$new(random_seed = -1)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{random_seed}}{(Optional) random seed for sampling}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{CppRNG} object.
}
}
}

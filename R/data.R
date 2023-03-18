#' Builtin primers
#'
#' Default forward and reverse primer of each 18S hypervariable regions, ITS1,
#' 5.8S, and ITS2.
#'
#' @format ## `builtin_primers`
#' A list with default primers:
#' \describe{
#'   \item{gene}{gene name, such as 18S, 5.8S, ITS}
#'   \item{region}{hypervariable region name, such as v1, v2, v3, v4, v5, v6,
#'   v7, v8, and v9 for 18S, 5.8S for 5.8S, ITS1 and ITS2 for ITS}
#'   \item{set}{only default set now}
#'   \item{name}{original name of the primer}
#'   \item{id}{primer id in pr2-primers database or other resources}
#'   \item{direction}{forward (fwd) or reverse (rev)}
#'   \item{specificity}{is the primer specific of a group}
#'   \item{seq}{primer sequence}
#'   \item{start_yeast}{start of primer relative to
#'   FU970071<https://www.ncbi.nlm.nih.gov/nuccore/FU970071>}
#'   \item{reference}{original reference where primer was first defined}
#' }
"builtin_primers"

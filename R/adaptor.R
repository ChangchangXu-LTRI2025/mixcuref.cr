#' Adapter for meta-layer: fit mixcuref.cr from a mc_spec
#'
#' The mixcuref.cr package's core functions work with a single survival formula
#' of the form \code{Surv(time, status) ~ x1 + x2 + ...} and assume the same
#' covariates appear in both the incidence and latency parts of the model.
#'
#' The meta-layer (e.g., Mixcuref.Meta) stores incidence and latency formulas
#' separately in a list-like \code{spec} object. This adapter builds a combined
#' survival formula using the union of covariate terms, constructs a default
#' initial value vector when not supplied, and dispatches to
#' \code{mixcure.penal.polar.plcr()}.
#'
#' @param spec A list-like object (typically class \code{mc_spec}) with at least
#'   \code{$data}, \code{$time}, \code{$status}, \code{$incidence}, and
#'   \code{$latency}.
#' @param control A list of engine-specific controls. Supported entries include:
#'   \itemize{
#'     \item \code{pl}: logical; use Firth-type penalized likelihood (default FALSE).
#'     \item \code{apct}: numeric; polar step size passed to \code{mixcure.penal.polar.plcr} (default 0.05).
#'     \item \code{init}: numeric vector of initial values (default constructed automatically).
#'   }
#' @return The fitted object returned by \code{mixcure.penal.polar.plcr()}.
#' @export
mixcurefcr_fit <- function(spec, control = list()) {
  stopifnot(is.list(spec))
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("mixcurefcr_fit(): package 'survival' is required.", call. = FALSE)
  }
  if (is.null(spec$data) || !is.data.frame(spec$data)) {
    stop("mixcurefcr_fit(): 'spec$data' must be a data.frame.", call. = FALSE)
  }

  dat <- spec$data
  time_var <- spec$time %||% "time"
  status_var <- spec$status %||% "status"

  if (!time_var %in% names(dat)) {
    stop("mixcurefcr_fit(): time column not found in data: ", time_var, call. = FALSE)
  }
  if (!status_var %in% names(dat)) {
    stop("mixcurefcr_fit(): status column not found in data: ", status_var, call. = FALSE)
  }

  # Build combined RHS from incidence + latency formulas.
  inc_terms <- character(0)
  lat_terms <- character(0)
  if (!is.null(spec$incidence) && inherits(spec$incidence, "formula")) {
    inc_terms <- attr(stats::terms(spec$incidence), "term.labels")
  }
  if (!is.null(spec$latency) && inherits(spec$latency, "formula")) {
    lat_terms <- attr(stats::terms(spec$latency), "term.labels")
  }
  rhs_terms <- unique(c(inc_terms, lat_terms))
  rhs <- if (length(rhs_terms) == 0) "1" else paste(rhs_terms, collapse = " + ")

  # Use survival::Surv(...) so we don't require attaching survival.
  formula <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", status_var, ") ~ ", rhs)
  )

  pl <- isTRUE(control$pl)
  apct <- control$apct %||% 0.05

  init <- control$init
  if (is.null(init)) {
    mf <- stats::model.frame(formula, data = dat, na.action = stats::na.omit)
    mm <- stats::model.matrix(formula, data = mf)
    k <- ncol(mm)
    init <- c(rep(0, 2 * k), 1)
  }

  mixcure.penal.polar.plcr(
    formula = formula,
    data = dat,
    apct = apct,
    init = init,
    pl = pl
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

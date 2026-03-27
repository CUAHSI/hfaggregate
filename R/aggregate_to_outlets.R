#' Partition flowpaths between outlets and flag mainstem by outlet levelpath
#'
#' @param flowpaths      data.frame or sf with flow network columns.
#' @param outlets        vector of outlet flowpath IDs (same type as id_col).
#' @param id_col         name of ID column (default "flowpath_id").
#' @param toid_col       name of downstream ID column (default "flowpath_toid").
#' @param levelpath_col  name of levelpath column (default "levelpathid").
#' @param hydroseq_col  name of hydroseq column (default "hydroseq").
#' @param include_outlet logical; if TRUE, the outlet feature belongs to its own group.
#' @param hydroseq_increases_downstream logical; if TRUE, higher hydroseq values are downstream.
#' @param auto_seed_terminals logical; if TRUE, terminal flowpaths are treated as outlets if no matches found.
#' @param quiet logical; if FALSE, prints summary messages.
#'
#' @return `flowpaths` with two new columns:
#'   - group_id (integer; NA where no downstream outlet is reachable)
#'   - mainstem (logical; TRUE to keep, FALSE to drop tributaries)

.partition_and_flag_mainstem <- function(flowpaths,
                                         outlets,
                                         id_col = "flowpath_id",
                                         toid_col = "flowpath_toid",
                                         levelpath_col = "levelpathid",
                                         hydroseq_col = "hydroseq",
                                         include_outlet = TRUE,
                                         hydroseq_increases_downstream = TRUE,
                                         auto_seed_terminals = FALSE,
                                         quiet = FALSE) {
  
  stopifnot(all(c(id_col, toid_col, levelpath_col, hydroseq_col) %in% names(flowpaths)))
  
  # --- Normalize & sanitize ----------------------------------------------------
  norm_chr <- function(x) trimws(as.character(x))
  ids      <- norm_chr(flowpaths[[id_col]])
  toids    <- norm_chr(flowpaths[[toid_col]])
  lvlpath  <- flowpaths[[levelpath_col]]
  hseq     <- flowpaths[[hydroseq_col]]
  outlets  <- norm_chr(outlets)
  
  # Treat common terminal codes as NA
  toids[toids %in% c("", "0", "NA", "NaN", "NULL", "None")] <- NA_character_
  
  # Build downstream index (NA if terminal or not found in ids)
  down_idx <- match(toids, ids)
  
  # Which rows are outlets (by id match)
  is_outlet <- ids %in% outlets
  n_match <- sum(is_outlet, na.rm = TRUE)
  
  if (n_match == 0L && !auto_seed_terminals) {
    bad <- setdiff(unique(outlets), unique(ids))
    stop(sprintf(
      "No provided outlets matched '%s'. Examples not found: %s",
      id_col,
      paste(utils::head(bad, 5), collapse = ", ")
    ))
  }
  
  # --- Seed labels at outlets (and optionally at terminals) --------------------
  nearest_out <- rep(NA_character_, length(ids))
  
  if (include_outlet && n_match > 0L) {
    nearest_out[is_outlet] <- ids[is_outlet]
  }
  if (auto_seed_terminals) {
    at_terminal <- is.na(down_idx)
    # Only seed terminals that don't already have an outlet label
    nearest_out[at_terminal & is.na(nearest_out)] <- ids[at_terminal & is.na(nearest_out)]
  }
  
  # If still no seeds at all, bail with actionable message
  if (all(is.na(nearest_out))) {
    stop("No seeds available to propagate (no outlet matches and auto_seed_terminals=FALSE).")
  }
  
  # --- Order-agnostic downstream→upstream propagation --------------------------
  changed <- TRUE; iters <- 0L
  while (changed) {
    iters <- iters + 1L
    prev <- nearest_out
    need <- is.na(nearest_out)                # rows that do not yet have a label
    di   <- down_idx[need]                    # their downstream index
    take <- !is.na(di)                        # only those with a valid downstream neighbor
    if (any(take)) {
      # Copy the downstream neighbor's current label
      nearest_out[which(need)[take]] <- nearest_out[di[take]]
    }
    changed <- !identical(prev, nearest_out)
    if (iters > 10000L) stop("Propagation did not converge; check for cycles in toids.")
  }
  
  # --- Build group_id keyed by outlet_id --------------------------------------
  uniq_outlets <- unique(nearest_out[!is.na(nearest_out)])
  group_id <- rep(NA_integer_, length(ids))
  if (length(uniq_outlets)) {
    group_id[!is.na(nearest_out)] <- match(nearest_out[!is.na(nearest_out)], uniq_outlets)
  }
  
  # --- Per-row outlet metadata -------------------------------------------------
  oi <- match(nearest_out, ids)                  # outlet row index for each feature
  outlet_levelpath <- ifelse(is.na(oi), NA, lvlpath[oi])
  outlet_hseq      <- ifelse(is.na(oi), NA, hseq[oi])
  
  # --- Mainstem: outlet levelpath *and* upstream-to-outlet by hydroseq --------
  if (isTRUE(hydroseq_increases_downstream)) {
    mainstem <- !is.na(group_id) & !is.na(outlet_levelpath) &
      (lvlpath == outlet_levelpath) &
      !is.na(hseq) & !is.na(outlet_hseq) &
      (hseq <= outlet_hseq)
  } else {
    mainstem <- !is.na(group_id) & !is.na(outlet_levelpath) &
      (lvlpath == outlet_levelpath) &
      !is.na(hseq) & !is.na(outlet_hseq) &
      (hseq >= outlet_hseq)
  }
  
  # --- Group hydroseq (copy the outlet's) -------------------------------------
  group_hydroseq <- outlet_hseq
  
  # --- Collapsed graph: effective downstream group ----------------------------
  outlet_row_by_group <- match(uniq_outlets, ids)
  outlet_down_idx     <- down_idx[outlet_row_by_group]
  next_outlet_id      <- ifelse(is.na(outlet_down_idx), NA_character_, nearest_out[outlet_down_idx])
  
  outlet_to_group <- setNames(seq_along(uniq_outlets), uniq_outlets)
  next_group_id   <- unname(outlet_to_group[next_outlet_id])
  
  effective_outlet_id <- ifelse(is.na(group_id), NA_character_, next_outlet_id[group_id])
  effective_toid      <- ifelse(is.na(group_id), NA_integer_,  next_group_id[group_id])
  
  # --- Optional debug summary --------------------------------------------------
  if (!quiet) {
    if (n_match == 0L && auto_seed_terminals) {
      message("Note: No provided outlets matched; groups were seeded from terminal flowpaths.")
    } else if (n_match > 0L) {
      message(sprintf("Seeded %d outlet(s); total groups formed: %d.", n_match, length(uniq_outlets)))
    }
  }
  
  # Bind back
  flowpaths$outlet_id            <- nearest_out
  flowpaths$group_id             <- group_id
  flowpaths$mainstem             <- mainstem
  flowpaths$group_hydroseq       <- group_hydroseq
  flowpaths$effective_outlet_id  <- effective_outlet_id
  flowpaths$effective_toid       <- effective_toid
  
  flowpaths
}


#' Aggregate flowpaths and divides to outlet groups
#'
#' Builds outlet-based groups, unions divides by group, and unions mainstem
#' flowpaths by group to produce an outlet-scale network layer. This function
#' identifies critical junctions along levelpaths carrying points-of-interest
#' (POIs) and treats their immediate upstream siblings as additional outlets,
#' ensuring tributary groups are preserved where POIs branch.
#'
#' @param gpkg Path to a hydrofabric GeoPackage.
#' @param vpu Optional VPU code to subset (if supported by your readers). Not
#'   currently used internally but kept for API compatibility.
#' @param flowpath Optional layer name for flowpaths in `gpkg`. If `NULL`,
#'   defaults used by `read_hydrofabric()`.
#' @param divide Optional layer name for divides in `gpkg`. If `NULL`, defaults
#'   used by `read_hydrofabric()`.
#' @param crs Target CRS used when reading the hydrofabric. Default: `5070`.
#' @param pois An `sf` (or data frame) with a `flowpath_id` column indicating
#'   POI-associated flowpaths to anchor mainstem levelpaths.
#' @param outlets Optional character/numeric vector of flowpath IDs to treat as
#'   outlets when `pois` is `NULL` or when additional outlets should be seeded.
#' @param outfile Optional output path (reserved for future write support).
#'
#' @return A named list with two `sf` layers:
#' \describe{
#'   \item{divides}{Polygons unioned by `group_id` (one per outlet group).}
#'   \item{flowpaths}{Mainstem lines unioned by `group_id` with attributes:
#'        `levelpathid`, `vpuid`, `member_comid` (CSV of members),
#'        and `flowpath_toid` representing the collapsed next-group link.}
#' }

#' @export
#' @importFrom dplyr filter select mutate group_by 
#' @importFrom dplyr left_join pull slice_max ungroup distinct rename
#' @importFrom sf st_drop_geometry st_as_sf

aggregate_to_outlets = function(gpkg = NULL,
                                vpu = NULL,
                                flowpath = NULL,
                                divide = NULL,
                                crs = 5070,
                                pois = NULL,
                                outlets = NULL,
                                outfile = NULL){
  
  id_col = "flowpath_id"
  toid_col = "flowpath_toid"
  levelpath_col = "levelpathid"
  hydroseq_col = "hydroseq"
  
  network_list <- read_hydrofabric(gpkg = gpkg, divides = divide,  flowpaths = flowpath, crs = crs) |>
    prepare_network() |>
    add_network_type()
  
  if (!is.null(pois)) {
    names(pois) <- tolower(names(pois))
    if (!"flowpath_id" %in% names(pois)) {
      cli::cli_abort("`pois` must contain a `flowpath_id` column.")
    }
    poi_outlets <- pois$flowpath_id
  } else {
    poi_outlets <- NULL
  }
  
  outlets <- unique(c(outlets, poi_outlets))
  outlets <- outlets[!is.na(outlets)]
  
  if (length(outlets) == 0) {
    cli::cli_abort("Provide either `pois` with `flowpath_id` or an `outlets` vector.")
  }
  
  mainstems_w_poi = filter(network_list$flowpaths, flowpath_id %in% outlets) |> 
    pull(levelpathid) |> 
    unique()
  
  critical_juntions <-  dplyr::select(network_list$flowpaths, 
                                      flowpath_toid = flowpath_id, 
                                      levelpathid_toid = levelpathid,
                                      hydroseq_toid = hydroseq) |> 
    sf::st_drop_geometry() |> 
    dplyr::distinct() |>
    dplyr::left_join(st_drop_geometry(network_list$flowpaths),
                     by = "flowpath_toid") |> 
    dplyr::select(flowpath_id, levelpathid, levelpathid_toid, hydroseq_toid) |>
    dplyr::group_by(flowpath_id) |> 
    dplyr::slice_max(hydroseq_toid) |> 
    dplyr::ungroup() |> 
    dplyr::filter(levelpathid != levelpathid_toid) |> 
    dplyr::right_join(dplyr::select(network_list$flowpaths, 
                                    flowpath_id, 
                                    flowpath_toid), 
                      by = "flowpath_id") |>
    sf::st_as_sf() |> 
    dplyr::filter(levelpathid %in% mainstems_w_poi) |> 
    node_geometry() |> 
    dplyr::rowwise() |>
    dplyr::mutate(
      ups_id = list(
        dplyr::filter(network_list$flowpaths,
                      flowpath_toid == .data$flowpath_toid,
                      flowpath_id != .data$flowpath_id) |>
          dplyr::pull(flowpath_id)
      )
    ) |>
    dplyr::ungroup()
  
  all_outlets  <- unique(c(unlist(critical_juntions$ups_id), outlets))
  critical_lps <- dplyr::pull(dplyr::filter(network_list$flowpaths, 
                                            flowpath_id %in% outlets), "levelpathid")
  
  fps          <- .partition_and_flag_mainstem(network_list$flowpaths, 
                                               all_outlets)
  
  divides_out <- network_list$divides |> 
    dplyr::select(flowpath_id) |> 
    dplyr::left_join(st_drop_geometry(select(fps, flowpath_id, group_id)), by = "flowpath_id") |>
    union_polygons("group_id") 
  
  tmp <-  dplyr::select(fps, group_id, member_comid, effective_toid) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(group_id) |>
    dplyr::mutate(member_comid = paste(member_comid, collapse = ","),
                  effective_toid = effective_toid[1]) |>
    dplyr::distinct() |>
    dplyr::ungroup()
  
  lps <- dplyr::filter(fps, mainstem) |>
    dplyr::distinct(group_id, levelpathid, vpuid)
  
  flowpaths_out <-  dplyr::filter(fps, levelpathid %in% critical_lps) |>
    union_linestrings("group_id") |>
    dplyr::left_join(lps, by = "group_id") |>
    dplyr::left_join(tmp, by = "group_id") |>
    dplyr::rename(flowpath_id = group_id, flowpath_toid = effective_toid)
  
  network_list = list(divides = divides_out, flowpaths = flowpaths_out)
  
  if(!is.null(pois)){
    network_list$pois = pois
  }

  
  if (!is.null(outfile)) {
    outfile <- write_hydrofabric(network_list,
                                 outfile,
                                 enforce_dm = FALSE)
    
    return(outfile)
  } else {
    return(list(divides = divides_out, flowpaths = flowpaths_out))
  }
}



#' Split a Vector at Given Positions
#'
#' @param x Vector.
#' @param pos Integer vector of split positions.
#' @return List of split vectors.
#' @noRd
splitAt <- function(x, pos) {
  unname(split(x, cumsum(seq_along(x) %in% (pos + 1))))
}
#' Assign Group IDs Based on Ideal Size
#'
#' Groups sequential values until their sum exceeds the threshold.
#' @param x Numeric vector.
#' @param athres Numeric. Threshold cumulative sum before starting new group.
#' @return Integer vector of group IDs.
#' @export
assign_id <- function(x, athres) {
  cumsum <- 0
  group <- 1
  result <- numeric()
  for (i in seq_along(x)) {
    cumsum <- cumsum + x[i]
    if (cumsum > athres) {
      group <- group + 1
      cumsum <- x[i]
    }
    result[i] <- group
  }
  return(result)
}

#' Merge Small Groups on the Edges
#'
#' Merges the first and last group if they're smaller than the threshold.
#' @param x Vector of values.
#' @param ind Grouping index for x.
#' @param thres Numeric threshold for merging.
#' @return Updated group index.
#' @export
pinch_sides <- function(x, ind, thres) {
  group_sums <- tapply(x, ind, sum)
  group_ids <- as.numeric(names(group_sums))
  if (length(group_sums) == 1) {
    return(ind)
  }
  if (group_sums[1] < thres) names(group_sums)[1] <- names(group_sums)[2]
  if (group_sums[length(group_sums)] < thres) names(group_sums)[length(group_sums)] <- names(group_sums)[length(group_sums) - 1]
  as.numeric(names(group_sums))[match(ind, group_ids)]
}

#' Merge Interior Groups Under a Threshold
#'
#' Merges internal groups that fall below a threshold with the smaller of their neighbors.
#' @param x Vector of values (area or length).
#' @param index_values Grouping index.
#' @param threshold Numeric. Minimum group sum threshold.
#' @return Updated group index.
#' @export
middle_massage <- function(x, index_values, threshold) {
  group_sums <- tapply(x, index_values, sum)
  group_ids <- as.numeric(names(group_sums))
  if (length(group_sums) == 1) {
    return(index_values)
  }
  small_groups <- which(group_sums < threshold)
  for (grp in small_groups) {
    if (grp == 1 || grp == length(group_sums)) next
    neighbors <- c(grp - 1, grp + 1)
    target <- neighbors[which.min(group_sums[neighbors])]
    names(group_sums)[grp] <- names(group_sums)[target]
  }
  as.numeric(names(group_sums))[match(index_values, group_ids)]
}

#' @title Aggregate Network to Uniform Size
#' @description This function aggregates a network to a desired size distribution while
#' enforcing minimum flowpath legnths and catchment areas. Additionally a set of explicit nexus
#' locations can be provided over which the network cannot be aggregated (see poi_to_outlet)
#' @param gpkg a path to a gpkg
#' @param divide If gpkg is NULL, then an sf data.frame, otherwise a the layer name. See details.
#' @param flowpath If gpkg is NULL, then an sf data.frame, otherwise a the layer name. See details.
#' @param outlets data.frame with mandatory "ID" column and optional "POI_ID" column. "ID" must be identifiers from
#' flowpath and divide data.frames and POI ID must be unique.
#' @param ideal_size_sqkm The ideal size of catchments (default = 10 sqkm)
#' @param min_length_km The minimum allowable length of flowpath features (default = 1 km)
#' @param min_area_sqkm The minimum allowable area of catchment features (default = 3 sqkm)
#' @param outfile of not NULL, where to write the output files
#' @param overwrite overwrite existing gf file. Default is FALSE
#' @param nexus_locations a data.frame with columns specifying the ID, and the nexus type.
#' @param log a filepath to write messages to or Boolean (TRUE = print to console; FALSE = no messages)
#' @param verbose print status updates. Default = TRUE
#' @return if outfile = TRUE, a file path, else a list object
#' @details If gpkg is not NULL, divide and flowpath can be left NULL as well. The code attempts to
#' infer the correct layers. The divides layer will be the one including the word "divide" or "catchment" and the
#' flowpath layer will be the one including 'flowpath' or 'flowline'. If no layers, or more then one layer are deemed possible
#' for each input, then the function will stop and ask for explicit names.
#' @export
#' @importFrom sf st_transform read_sf st_set_crs write_sf st_layers st_is_valid
#' @importFrom dplyr left_join filter semi_join
#' @importFrom tidyr separate_longer_delim
#' @importFrom hfutils clean_geometry

aggregate_to_distribution <- function(gpkg = NULL,
                                      vpu = NULL,
                                      flowpath = NULL,
                                      divide = NULL,
                                      crs = 5070,
                                      pois = NULL,
                                      max_junction = 2,
                                      ideal_size_sqkm = 10,
                                      min_length_km = 1,
                                      min_area_sqkm = 3,
                                      outfile = NULL,
                                      overwrite = FALSE,
                                      verbose = TRUE) {

  if (!is.null(outfile)) {
    if (file.exists(outfile) & overwrite) {
      unlink(outfile)
    } else if (file.exists(outfile)) {
      cli::cli_alert_warning(glue("{outfile} already exists and overwrite is FALSE"))
      return(outfile)
    }
  }

  network_list <- read_hydrofabric(gpkg = gpkg,
                                   divides = divide, 
                                   flowpaths = flowpath,
                                   crs = crs) |>
    prepare_network() |>
    add_network_type(verbose = FALSE)
  
  if (!"member_comid" %in% names(network_list$flowpaths)) {
    network_list$flowpaths$member_comid <- network_list$flowpaths$flowpath_id
  }
  
  has_pois <- !is.null(pois)
  
  warn_poi <- function(stage, values) {
    if (has_pois) {
      total <- length(unique(pois$poi_id))
      tracked <- length(unique(na.omit(values)))
      if (!identical(total, tracked)) {
        warning(sprintf("All POIs are not indexable: %s", stage), call. = FALSE)
      }
    }
  }
  
  # Add outlets
  if (has_pois) {
    names(pois) <- tolower(names(pois))
    
    if (!all(c("flowpath_id", "poi_id") %in% names(pois))) {
      cli::cli_abort("`pois` must include `flowpath_id` and `poi_id` columns.")
    }
    
    bad <- filter(pois, !flowpath_id %in% network_list$flowpaths$flowpath_id)
    
    outflows <- pois %>%
      st_drop_geometry() %>%
      select(poi_id, flowpath_id)
    
    network_list$flowpaths <- left_join(mutate(network_list$flowpaths, poi_id = NULL),
                                        st_drop_geometry(outflows),
                                        by = "flowpath_id"
    )
    
    warn_poi("reference", network_list$flowpaths$poi_id)
  } else {
    network_list$flowpaths$poi_id <- NA
    network_list$flowpaths$hl_uri <- NA
    outflows <- NULL
  }
  
  main_agg <- aggregate_along_mainstems(
    network_list,
    max_junction = max_junction,
    ideal_size_sqkm,
    min_area_sqkm,
    min_length_km)
  
  warn_poi("mainstem aggregation", main_agg$flowpaths$poi_id)
  
  collapse_agg <- collapse_headwaters(network_list = main_agg, min_area_sqkm, min_length_km)
  
  collapse_agg$divides <- clean_geometry(collapse_agg$divides,  ID = "divide_id", keep = NULL)
  
  warn_poi("collapse", collapse_agg$flowpaths$poi_id)
  
  network_list <- collapse_agg
  rm(collapse_agg)
  rm(main_agg)
  
  if (!is.null(pois)) {
    network_list$pois <- pois
  }
  
  network_list$flowpaths <-
    select(
      network_list$flowpaths,
      flowpath_id,
      flowpath_toid,
      flowline_id,
      mainstem = levelpathid,
      member_comid,
      poi_id,
      hydroseq,
      lengthkm,
      areasqkm,
      tot_drainage_areasqkm = tot_drainage_area,
      has_divide
    ) %>%
    mutate(divide_id = ifelse(flowpath_id %in% network_list$divides$divide_id, flowpath_id, NA))
  
  warn_poi("final", network_list$flowpaths$poi_id)
  
  topo <- st_drop_geometry(network_list$flowpaths) %>%
    select(divide_id, flowpath_toid)
  
  network_list$divides <- select(network_list$divides, divide_id, areasqkm) %>%
    mutate(
      divide_id,
      has_flowline = TRUE,
      ds_id = NA,
      type = "network"
    ) %>%
    left_join(topo, by = "divide_id")
  
  network_list$network <- st_drop_geometry(network_list$flowpaths) %>%
    select(
      flowpath_id,
      flowpath_toid,
      member = member_comid,
      flowline_id,
      divide_id,
      any_of("poi_id"),
      mainstem,
      hydroseq,
      lengthkm,
      areasqkm,
      tot_drainage_areasqkm) %>%
    separate_longer_delim(col = "member", delim = ",") %>%
    mutate(
      hf_id_part = sapply(
        strsplit(member, "[.]"),
        FUN = function(x) {
          x[2]
        }
      ),
      hf_id_part = ifelse(is.na(hf_id_part), 1L, floor(as.numeric((
        hf_id_part
      )))),
      hf_id = sapply(
        strsplit(member, "[.]"),
        FUN = function(x) {
          as.numeric(x[1])
        }
      ),
      member = NULL,
      hf_source = "NHDPlusV2"
    ) %>%
    left_join(st_drop_geometry(select(
      network_list$divides, divide_id, type, ds_id
    )), by = "divide_id") |>
    distinct()
  
  warn_poi("network", network_list$network$poi_id)
  
  if (!is.null(vpu)) {
    network_list$network$vpu <- vpu
  } else {
    network_list$network$vpu <- NA
  }
  
  if (!all(st_geometry_type(network_list$divides) == "POLYGON")) {
    warning("MULTIPOLYGONS FOUND VPU: ", vpu)
  }
  
  if (!is.null(outfile)) {
    outfile <- write_hydrofabric(network_list,
                                 outfile,
                                 verbose = verbose,
                                 enforce_dm = FALSE
    )
    
    return(outfile)
  } else {
    network_list
  }
}

#' Aggregate Catchments Along Mainstems
#'
#' Groups small or short flowpaths and their associated divides into larger units
#' based on ideal area and length thresholds. Ensures topology is preserved.
#' @param network_list A list of `flowpaths` and `divides` as sf objects.
#' @param ideal_size_sqkm Numeric. Target size for catchments.
#' @param min_area_sqkm Numeric. Minimum catchment area.
#' @param min_length_km Numeric. Minimum flowpath length.
#' @param verbose Logical. Whether to print log messages.
#' @return Updated list of flowpaths and divides aggregated by size and topology.
#' @importFrom glue glue
#' @importFrom dplyr select distinct filter group_by slice ungroup mutate summarise cur_group_id n arrange slice_min

aggregate_along_mainstems <- function(network_list,
                                      max_junction,
                                      ideal_size_sqkm,
                                      min_area_sqkm,
                                      min_length_km) {
  
  cli::cli_alert_info("\n---  Aggregate Along Mainstem ---")
  cli::cli_alert_info(glue("max junction --> {max_junction}"))
  cli::cli_alert_info(glue("ideal_size_sqkm --> {ideal_size_sqkm}"))
  cli::cli_alert_info(glue("min_length_km --> {min_length_km}"))
  cli::cli_alert_info(glue("min_area_sqkm --> {min_area_sqkm}"))
  
  upstream_index <- network_list$flowpaths |>
    st_drop_geometry() |>
    select(flowpath_id, hl_un = poi_id) |>
    distinct() |>
    filter(!is.na(hl_un)) |>
    group_by(flowpath_id) |>
    slice(1) |>
    ungroup()
  
  fline <- network_list$flowpaths |>
    st_drop_geometry() |>
    mutate(hl_dn = poi_id) |>
    left_join(upstream_index, by = "flowpath_id") |>
    distinct()
  
  index_table <- fline |>
    group_by(levelpathid) |>
    mutate(ind = cs_group(areasqkm, lengthkm, max_junction, hl_dn, hl_un, ideal_size_sqkm, min_area_sqkm, min_length_km)) |>
    ungroup() |>
    group_by(levelpathid, ind) |>
    mutate(set = cur_group_id(), n = n()) |>
    ungroup() |>
    select(set, flowpath_id, flowpath_toid, levelpathid, hydroseq, member_comid, poi_id, n)
  
  aggregated <- aggregate_sets(network_list, index_table)
  aggregated <- add_network_type(aggregated, verbose = verbose)
  
  aggregated$lookup <- index_table |>
    group_by(set) |>
    summarise(flowline_id = paste(flowpath_id, collapse = ","), 
              flowline_toid = paste(flowpath_toid, collapse = ","),
              member_comid = paste(member_comid, collapse = ","), .groups = "drop") |>
    rename(flowpath_id = set)
  
  cli::cli_alert_success(glue("Merged to idealized catchment size of {ideal_size_sqkm} sqkm: {nrow(network_list$flowpaths) - nrow(aggregated$flowpaths)} features removed"))

  return(aggregated)
}

#' Aggregate Sets by Index Table
#' @param network_list a list of flowpaths and catchments
#' @param index_table index table to aggregate with
#' @return a list of catchments and flowpaths that have been validated
#' @importFrom dplyr group_by mutate slice_max ungroup select left_join everything filter bind_rows rename %>% inner_join
#' @importFrom sf st_as_sf
#' @importFrom hfutils rename_geometry union_polygons union_linestrings add_areasqkm write_hydrofabric
#' @importFrom terra vect combineGeoms

aggregate_sets = function(network_list, index_table) {
  
  set_topo = index_table %>%
    group_by(set) %>%
    mutate(
      member_comid  = paste(member_comid, collapse = ","),
      poi_id  = paste(poi_id[!is.na(poi_id)], collapse = ","),
      poi_id  = ifelse(poi_id == "", NA, poi_id)
    ) %>%
    arrange(hydroseq) %>%
    select(set, flowpath_id, flowpath_toid, levelpathid, poi_id, hydroseq, member_comid) %>%
    ungroup()
  
  set_topo_fin = left_join(
    select(
      set_topo,
      set,
      flowpath_id = flowpath_toid,
      hydroseq,
      levelpathid,
      poi_id,
      member_comid
    ),
    select(set_topo, toset = set, flowpath_id),
    by = "flowpath_id"
  ) %>%
    group_by(set) %>%
    mutate(toset = ifelse(is.na(toset), 0, toset)) %>%
    filter(set != toset) %>%
    slice_min(hydroseq) %>%
    ungroup() %>%
    select(set, toset, levelpathid, poi_id, member_comid)
  
  ####
  
  single_flowpaths = tryCatch({
    filter(index_table, n == 1) %>%
      inner_join(select(network_list$flowpaths, flowpath_id), by = "flowpath_id") %>%
      st_as_sf() %>%
      select(set) %>%
      rename_geometry("geometry")
  }, error = function(e) {
    NULL
  })
  
  flowpaths_out  = tryCatch({
    filter(index_table, n > 1) %>%
      inner_join(network_list$flowpaths, by = "flowpath_id") %>%
      st_as_sf() %>%
      select(set) %>%
      union_linestrings('set') %>%
      rename_geometry("geometry")
  }, error = function(e) {
    NULL
  })
  
  flowpaths_out = flowpaths_out |>
    bind_rows(single_flowpaths) %>%
    left_join(set_topo_fin, by = "set") %>%
    rename(flowpath_id = set, flowpath_toid = toset)
  
  ####
  single_catchments = tryCatch({
    filter(index_table, n == 1) %>%
      inner_join(network_list$divides, by = "flowpath_id") %>%
      st_as_sf() %>%
      select(set) %>%
      rename_geometry("geometry")
  }, error = function(e) {
    NULL
  })
  
  catchments_out  = tryCatch({
    filter(index_table, n != 1) %>%
      inner_join(network_list$divides, by = "flowpath_id") %>%
      st_as_sf() %>%
      select(set) %>%
      union_polygons('set') %>%
      mutate(areasqkm = add_areasqkm(.))
  }, error = function(e) {
    NULL
  })
  
  catchments_out = catchments_out %>%
    bind_rows(single_catchments) %>%
    left_join(set_topo_fin, by = "set") %>%
    select(divide_id = set, flowpath_toid = toset)
  
  ####
  
  if (!is.null(catchments_out)) {
    
    mps = suppressWarnings({
      catchments_out %>%
        st_cast("MULTIPOLYGON") %>%
        st_cast("POLYGON") %>%
        sf::st_make_valid() %>%
        mutate(n = n(), by = divide_id)
    })
    
    if (nrow(mps) > nrow(catchments_out)) {
      fixers = filter(mps, n > 1) %>%
        mutate(areasqkm = add_areasqkm(.), tmpID = 1:n()) %>%
        group_by(divide_id)
      
      keep = slice_max(fixers, areasqkm) %>%
        ungroup()
      
      to_merge = filter(fixers, !tmpID %in% keep$tmpID) %>%
        ungroup()
      
      good_to_go = filter(mps, n == 1) %>%
        bind_rows(select(keep, divide_id, flowpath_toid)) |> 
        sf::st_make_valid()
      

      catchments_out = suppressWarnings({
        terra::combineGeoms(terra::vect(good_to_go), terra::vect(to_merge)) %>%
          st_as_sf() %>%
          st_cast("POLYGON")
      })
    }
  }
  
  prepare_network(network_list = list(flowpaths = flowpaths_out, divides = catchments_out))
}

#' Group Flowpaths by Area and Length Constraints
#'
#' Groups flowpaths along a levelpath by cumulative area and length thresholds,
#' accounting for enforced breaks at upstream/downstream POIs.
#' @param areas Numeric vector of flowpath areas (sqkm).
#' @param lengths Numeric vector of flowpath lengths (km).
#' @param exclude_dn Vector of downstream break indicators.
#' @param exclude_un Vector of upstream break indicators.
#' @param ideal_size_sqkm Numeric. Target group area size.
#' @param amin Minimum allowed group area.
#' @param lmin Minimum allowed group length.
#' @return Integer vector assigning group IDs.
#' @export

cs_group <- function(areas, lengths, max_junction, exclude_dn, exclude_un,
                     ideal_size_sqkm, amin, lmin) {
  
  areas[is.na(areas)]     <- 0
  
  lengths[is.na(lengths)] <- 0
  
  if (length(areas) == 1) return(1)
  
  # Build breaks from exclude flags
  breaks <- sort(c(which(!is.na(exclude_dn)), which(!is.na(exclude_un)) - 1))
  
  area_chunks   <- splitAt(areas,   breaks)
  length_chunks <- splitAt(lengths, breaks)
  stopifnot(length(area_chunks) == length(length_chunks))
  
  # Initial grouping by ideal area
  initial <- lapply(area_chunks, assign_id, athres = ideal_size_sqkm)
  
  # "Pinch" ends: meet area then length minima
  pinched_area <- mapply(pinch_sides,  area_chunks,   initial,
                         MoreArgs = list(thres = amin), SIMPLIFY = FALSE)
  pinched_len  <- mapply(pinch_sides,  length_chunks, pinched_area,
                         MoreArgs = list(thres = lmin), SIMPLIFY = FALSE)
  
  # "Massage" middles: meet area then length minima
  massaged_area <- mapply(middle_massage, area_chunks,   pinched_len,
                          MoreArgs = list(threshold = amin), SIMPLIFY = FALSE)
  massaged_len  <- mapply(middle_massage, length_chunks, massaged_area,
                          MoreArgs = list(threshold = lmin), SIMPLIFY = FALSE)
  
  # ---- NEW: enforce max_junction per chunk ----------------------------------
  enforce_max_junction <- function(ids, max_junction) {
    # No cap if max_junction is NA/NULL/Inf
    if (is.null(max_junction) || is.na(max_junction) || is.infinite(max_junction)) {
      return(ids)
    }
    # Minimum sensible: 0 junctions => groups of size 1
    max_features <- max(1, max_junction + 1)
    
    r <- rle(ids)
    new_ids <- integer(length(ids))
    idx <- 1L
    local_group_counter <- 0
    
    for (k in seq_along(r$lengths)) {
      L <- r$lengths[k]
      # Allowed max features per group = max_junction + 1
      nseg <- ceiling(L / max_features)
      for (s in seq_len(nseg)) {
        len_s <- if (s < nseg) max_features else L - max_features * (nseg - 1)
        local_group_counter <- local_group_counter + 1
        new_ids[idx:(idx + len_s - 1)] <- local_group_counter
        idx <- idx + len_s
      }
    }
    new_ids
  }
  
  massaged_len <- lapply(massaged_len, enforce_max_junction, max_junction = max_junction)
  # ----------------------------------------------------------------------------
  
  # Keep labels unique across chunks
  for (i in seq_along(massaged_len)) {
    massaged_len[[i]] <- massaged_len[[i]] + 1e9 * i
  }
  
  unlist(massaged_len)
}


collapse_headwaters <- function(network_list,
                                 min_area_sqkm = 3,
                                 min_length_km = 1) {
  
  cli::cli_alert_info("\n--- Collapse Network Inward ---\n")
  
  network_list$flowpaths <- network_list$flowpaths |>
    left_join(select(network_list$lookup, flowpath_id, flowline_id, flowline_toid), by = "flowpath_id")
  
  start <- nrow(network_list$flowpaths)
  
  mapping_table <- build_headwater_collapse(network_list, 
                                            min_area_sqkm, 
                                            min_length_km)
  
  count <- 0
  
  while (nrow(mapping_table) > 0) {
    count <- count + 1
    
    cli::cli_alert_info(glue("Collapsing: {nrow(mapping_table)} features (round {count})"), verbose)
    
    ind <- match(mapping_table$becomes, network_list$flowpaths$flowpath_id)
    
    network_list$flowpaths$member_comid[ind] <- mapping_table$member_comid
    network_list$flowpaths$flowline_id[ind] <- mapping_table$fl_id
    network_list$flowpaths$flowline_toid[ind] <- mapping_table$fl_id
    
    fl <- filter(network_list$flowpaths, !flowpath_id %in% mapping_table$flowpath_id)
    
    cat <- filter(
      network_list$divides,
      divide_id %in% c(mapping_table$flowpath_id, mapping_table$becomes)
    ) |>
      left_join(mapping_table, by = c("divide_id" = "flowpath_id")) |>
      mutate(divide_id = ifelse(is.na(becomes), divide_id, becomes)) |>
      union_polygons("divide_id") |>
      bind_rows(filter(
        select(network_list$divides, divide_id),
        !divide_id %in% c(mapping_table$flowpath_id, mapping_table$becomes)
      )) %>%
      left_join(st_drop_geometry(select(fl, flowpath_id, flowpath_toid)), by = c("divide_id" = "flowpath_id"))
    
    network_list <- prepare_network(network_list = list(flowpaths = fl, divides = cat))
    
    mapping_table <- build_headwater_collapse(network_list, min_area_sqkm, min_length_km)
  }
  
  network_list <- add_network_type(network_list, verbose)
  
  cli::cli_alert_success(glue("Collapsed {start - nrow(network_list$flowpaths)} features."))
  
  return(network_list)
}

build_headwater_collapse <- function(network_list, 
                                     min_area_sqkm = 3, 
                                     min_length_km = 1) {
  
  touch_id <- define_touch_id(flowpaths = network_list$flowpaths) %>%
    filter(flowpath_id != touches) %>%
    group_by(flowpath_id) %>%
    mutate(check = ifelse(flowpath_toid == touches, 1, 2)) %>%
    slice_min(check) %>%
    ungroup() %>%
    select(-flowpath_toid)
  
  divide_touch_id <- st_intersects(network_list$divides)
  
  divide_touches <- data.frame(
    divide_id = rep(network_list$divides$divide_id, times = lengths(divide_touch_id)),
    becomes = network_list$divides$divide_id[unlist(divide_touch_id)]
  )
  
  # bad fps are those that are both hw and too small or too large
  df <- left_join(network_list$flowpaths, select(touch_id, -poi_id), by = "flowpath_id") %>%
    mutate(
      inflow = ifelse(flowpath_id %in% touch_id$touches, TRUE, FALSE),
      hw = ifelse(!flowpath_id %in% flowpath_toid, TRUE, FALSE),
      hw = ifelse(hw & !inflow, TRUE, FALSE),
      small = areasqkm < min_area_sqkm | lengthkm < min_length_km
    ) %>%
    filter(hw, small) %>%
    st_drop_geometry() %>%
    select(flowpath_id, becomes = touches, member_comid, poi_id) %>%
    filter(becomes != 0)
  
  df$mC1 <- network_list$flowpaths$member_comid[match(df$flowpath_id, network_list$flowpaths$flowpath_id)]
  df$mC2 <- network_list$flowpaths$member_comid[match(df$becomes, network_list$flowpaths$flowpath_id)]
  
  df$hl1 <- network_list$flowpaths$poi_id[match(df$flowpath_id, network_list$flowpaths$flowpath_id)]
  df$hl2 <- network_list$flowpaths$poi_id[match(df$becomes, network_list$flowpaths$flowpath_id)]
  
  df$fl1 <- network_list$flowpaths$flowline_id[match(df$flowpath_id, network_list$flowpaths$flowpath_id)]
  df$fl2 <- network_list$flowpaths$flowline_id[match(df$becomes, network_list$flowpaths$flowpath_id)]
  
  tmp <- suppressWarnings({
    group_by(df, becomes) %>%
      mutate(
        member_comid = paste0(mC2[1], ",", paste(mC1, collapse = ",")),
        poi_id = as.numeric(paste(unique(na.omit(c(hl1, hl2))), collapse = ",")),
        fl_id = paste0(fl2[1], ",", paste(fl1, collapse = ","))
      ) %>%
      ungroup() %>%
      select(-mC1, -mC2, -hl1, -hl2, -fl1, -fl2)
  })
  
  
  semi_join(tmp, rename(divide_touches, flowpath_id = divide_id), by = c("flowpath_id", "becomes"))
}



collapse_headwaters <- function(network_list,
                                 min_area_sqkm = 3,
                                 min_length_km = 1,
                                 verbose = TRUE,
                                 cache_file = NULL) {
  
  cli::cli_alert_info("\n--- Collapse Network Inward ---\n")
  
  network_list$flowpaths <- network_list$flowpaths |>
    left_join(select(network_list$lookup, flowpath_id, flowline_id), by = "flowpath_id")
  
  start <- nrow(network_list$flowpaths)
  
  mapping_table <- build_headwater_collapse(network_list, 
                                            min_area_sqkm, 
                                            min_length_km)
  
  count <- 0
  
  while (nrow(mapping_table) > 0) {
    count <- count + 1
    
    cli::cli_alert_info(glue("Collapsing: {nrow(mapping_table)} features (round {count})"), verbose)
    
    ind <- match(mapping_table$becomes, network_list$flowpaths$flowpath_id)
    
    network_list$flowpaths$member_comid[ind] <- mapping_table$member_comid
    network_list$flowpaths$flowline_id[ind] <- mapping_table$fl_id
    
    fl <- filter(network_list$flowpaths, !flowpath_id %in% mapping_table$flowpath_id)
    
    cat <- filter(
      network_list$divides,
      divide_id %in% c(mapping_table$flowpath_id, mapping_table$becomes)
    ) |>
      left_join(mapping_table, by = c("divide_id" = "flowpath_id")) |>
      mutate(divide_id = ifelse(is.na(becomes), divide_id, becomes)) |>
      union_polygons("divide_id") |>
      bind_rows(filter(
        select(network_list$divides, divide_id),
        !divide_id %in% c(mapping_table$flowpath_id, mapping_table$becomes)
      )) %>%
      left_join(st_drop_geometry(select(fl, flowpath_id, flowpath_toid)), by = c("divide_id" = "flowpath_id"))
    
    network_list <- prepare_network(network_list = list(flowpaths = fl, divides = cat))
    
    mapping_table <- build_headwater_collapse(network_list, min_area_sqkm, min_length_km)
  }
  
  network_list <- add_network_type(network_list, verbose)
  
  cli::cli_alert_success(glue("Collapsed {start - nrow(network_list$flowpaths)} features."))
  
  return(network_list)
}


#' Identify intersection types and downstream topology
#' @param flowpaths sf LINESTRING
#' @return data.frame with id, type, touches, touches_toID columns
#' @export
#' @importFrom sf st_intersects st_drop_geometry st_cast
#' @importFrom dplyr left_join select
#' @importFrom hfutils get_node node_geometry

define_touch_id <- function(flowpaths, term_cut = 1e9) {
  
  tmp <- st_cast(flowpaths, "MULTILINESTRING")
  
  ends <- rename_geometry(tmp, "geometry") |> 
    mutate(geometry = get_node(tmp, "end"))
  
  starts_ends <- bind_rows(data.frame(flowpath_id = tmp$flowpath_id, get_node(tmp, "start")), 
                           data.frame(flowpath_id = tmp$flowpath_id, get_node(tmp, "end"))) |> 
    st_as_sf()
  
  emap <- st_intersects(ends, starts_ends)
  tmp$type <- ifelse(lengths(emap) > 1, "nex", "jun")
  tmp$type <- ifelse(tmp$flowpath_toid > term_cut, "term", tmp$type)
  
  ends2 <- left_join(st_drop_geometry(select(ends, flowpath_id)),
                     st_drop_geometry(select(tmp, flowpath_id, flowpath_toid, type, hydroseq, poi_id)),
                     by = "flowpath_id"
  )
  
  tmap <- st_intersects(ends, tmp)
  
  data.frame(
    flowpath_id = rep(ends2$flowpath_id, times = lengths(tmap)),
    flowpath_toid = rep(ends2$flowpath_toid, times = lengths(tmap)),
    type = rep(ends2$type, times = lengths(tmap)),
    hs = rep(ends2$hydroseq, times = lengths(tmap)),
    poi_id = rep(ends2$poi_id, times = lengths(tmap)),
    touches = tmp$flowpath_id[unlist(tmap)],
    touches_toID = tmp$flowpath_toid[unlist(tmap)],
    touches_hs = tmp$hydroseq[unlist(tmap)]
  ) %>%
    filter(is.na(poi_id))
}

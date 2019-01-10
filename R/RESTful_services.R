
# *************************************************************
# RESTful web service with data as JSON arrays ====
# *************************************************************

json2script = function(json) {
  script = json$code$source
  script = gsub("\"", "'", script)
  script = gsub("\r", "", script)
  # Todo: Source the string `script` directly without writing it to disk
  # Note: Need to remove escape characters in `script`

  print(Sys.time())
  cat("Extracted script for JSON!\n")
  
  return(script)
}

tile2raster = function(tile, time_num, proj) {
  xmin = tile$extent$west
  xmax = tile$extent$east
  ymin = tile$extent$south
  ymax = tile$extent$north
  resx = tile$extent$width
  resy = tile$extent$height

  xtot = length(tile$data[[1]])
  ytot = length(tile$data[[1]][[1]])

  xyz = matrix(ncol = 3, nrow = xtot * ytot)
  xyz = as.data.frame(xyz)
  colnames(xyz) = c("x", "y","z")

  x = xmin + resx/2
  y = ymin + resy/2
  xyz$x = rep(seq(from = x, to = xmax - resx/2, by = resx), each = ytot)
  xyz$y = rep(seq(from = y, to = ymax - resy/2, by = resy), xtot)
  xyz$z = rapply(object = tile$data[[time_num]], f = unlist)

  return(raster::rasterFromXYZ(xyz, crs = proj))
}

json2stars = function(json) {
  print(Sys.time())
  cat("Started converting JSON to stars...\n")
  proj_string = json$data$proj
  num_bands = length(json$data$raster_collection_tiles)
  num_time = length(json$data$raster_collection_tiles[[1]]$start_times)

  bt_list = list()
  length(bt_list) = num_time

  timestamps = strptime(json$data$raster_collection_tiles[[1]]$start_times, format = "%Y-%m-%dT%T", tz = "Europe/Berlin")
  
  # #Todo: Get rid of the nested for loops below and in `tile2raster()`
  # #Todo: Think about using dataframes/martices

  on_bands = function(band_num, raster_collection, time_num, proj_string) {
    cat(paste("\n", Sys.time(), "; Processing Band: ", band_num, "; ", sep = ""))
    tile = raster_collection[[band_num]]
    
    return(tile2raster(tile = tile, time_num = time_num, proj = proj_string))
  }

  on_times = function(time_num, raster_collection, proj_string) {
    cat(paste("\n", Sys.time(), "; Time: ", time_num, "...", sep = ""))
    
    return(lapply(
      X = as.list(1:num_bands),
      FUN = on_bands,
      raster_collection,
      time_num,
      proj_string
    ))
  }

  raster_collection = json$data$raster_collection_tiles
  bt_list = lapply(X = as.list(1:num_time), FUN = on_times, raster_collection, proj_string)
  cat("\n")
  print(Sys.time())
  cat("Finished converting JSON to Raster objects!\n\n")

  stars_obj = NULL
  print(Sys.time())
  cat("Starting to convert Raster to stars objects...\n")
  stars_bands = list()
  for(time_num in 1:num_time) {
    as_stars = lapply(X = bt_list[[time_num]], FUN = stars::st_as_stars)

    append_stars = append(as_stars, values = c(along = "band"))
    stars_bands[[time_num]] = do.call(c, append_stars)
  }
  stars_obj = append(stars_bands, values = c(along = "time"))
  stars_obj = do.call(c, stars_obj)
  attr(stars_obj, "dimensions")$time$values = timestamps
  attr(stars_obj, "dimensions")$time$offset = timestamps[1]
  attr(stars_obj, "dimensions")$time$delta = mean(diff(timestamps))


  print(Sys.time())
  cat("Converted JSON to stars object!\n\n")
  stars_obj
}

run_script = function(stars_obj, dim_mod, script_text) {
  print(Sys.time())
  cat("Started executing UDF on stars object...\n")
  # dim_mod = 1, 2 means space
  # dim_mod = 3 means band
  # dim_mod = 4 means time
  # dim_mod = 5 means whether raster or feature (default: raster)
  in_dim = dim(stars_obj)
  all_dim = 1:4
  if("x" %in% names(in_dim) && "y" %in% names(in_dim)) all_dim[1] = 1 else all_dim[1] = NA
  
  if("band" %in% names(in_dim)) all_dim[2] = 2 else all_dim[2] = NA
  
  if("time" %in% names(in_dim)) all_dim[3] = 3 else all_dim[3] = NA
  
  all_dim[4] = 4 #Currently assuming `stars_obj` has rasters
  parsed_script = parse(text = script_text)
  
  if(is.expression(parsed_script)) {
    function_name = eval(parsed_script)
    result = stars::st_apply(stars_obj, FUN = function_name, MARGIN = all_dim[-c(dim_mod)])
    new_dim = all_dim
    new_dim[dim_mod] = NA
  } else stop("Script text is unavailable or is not a valid expression!")
  
  print(Sys.time())
  cat("Applied UDF on stars object!\n\n")
  
  return(result)
}

run_script_raw = function(stars_obj, script_text) {
  parsed_script = parse(text = script_text)
  
  if (is.expression(parsed_script)) {
    function_name = eval(parsed_script)
    result = function_name(stars_obj)
  } else stop("Script text is unavailable or is not a valid expression!")
  
  if (class(result) == "stars") return(result) else return(stars_obj)
}

stars2json = function(stars_obj, json_in) {
  print(Sys.time())
  cat("Started converting stars object to JSON...\n\n")
  json_out = json_in[["data"]] #Copying structure of JSON but only the element "data"
  json_out$proj = attr(stars_obj, "dimensions")$x$refsys
  tot_bands = as.numeric(dim(stars_obj)["band"])

  calc_extent = function(stars_obj, bands) {
    if(!is.na(bands)) {
      # Need another robust way to loop over bands & time since using `attr()` in the manner
      # below will not work for stars objects with arbitrary dimensions
      if("time" %in% dimnames(stars_obj)) {
        delta_x = attr(stars_obj[,,,bands,], "dimensions")$x$delta
        delta_y = attr(stars_obj[,,,bands,], "dimensions")$y$delta
        x1 = attr(stars_obj[,,,bands,], "dimensions")$x$offset
        x2 = attr(stars_obj[,,,bands,], "dimensions")$x$offset + delta_x * attr(stars_obj[,,,bands,], "dimensions")$x$to
        y1 = attr(stars_obj[,,,bands,], "dimensions")$y$offset
        y2 = attr(stars_obj[,,,bands,], "dimensions")$y$offset + delta_y * attr(stars_obj[,,,bands,], "dimensions")$y$to
      } else {
        delta_x = attr(stars_obj[,,,bands], "dimensions")$x$delta
        delta_y = attr(stars_obj[,,,bands], "dimensions")$y$delta
        x1 = attr(stars_obj[,,,bands], "dimensions")$x$offset
        x2 = attr(stars_obj[,,,bands], "dimensions")$x$offset + delta_x * attr(stars_obj[,,,bands], "dimensions")$x$to
        y1 = attr(stars_obj[,,,bands], "dimensions")$y$offset
        y2 = attr(stars_obj[,,,bands], "dimensions")$y$offset + delta_y * attr(stars_obj[,,,bands], "dimensions")$y$to
      }
    } else {
      if("time" %in% dimnames(stars_obj)) {
        delta_x = attr(stars_obj[,,,], "dimensions")$x$delta
        delta_y = attr(stars_obj[,,,], "dimensions")$y$delta
        x1 = attr(stars_obj[,,,], "dimensions")$x$offset
        x2 = attr(stars_obj[,,,], "dimensions")$x$offset + delta_x * attr(stars_obj[,,,], "dimensions")$x$to
        y1 = attr(stars_obj[,,,], "dimensions")$y$offset
        y2 = attr(stars_obj[,,,], "dimensions")$y$offset + delta_y * attr(stars_obj[,,,], "dimensions")$y$to
      } else {
        delta_x = attr(stars_obj[,,], "dimensions")$x$delta
        delta_y = attr(stars_obj[,,], "dimensions")$y$delta
        x1 = attr(stars_obj[,,], "dimensions")$x$offset
        x2 = attr(stars_obj[,,], "dimensions")$x$offset + delta_x * attr(stars_obj[,,], "dimensions")$x$to
        y1 = attr(stars_obj[,,], "dimensions")$y$offset
        y2 = attr(stars_obj[,,], "dimensions")$y$offset + delta_y * attr(stars_obj[,,], "dimensions")$y$to
      }
    }
    
    return(
      list(north = max(y1,y2), 
           south = min(y1,y2), 
           west = min(x1,x2), 
           east = max(x1,x2), 
           height = if(sign(delta_y) < 0) -1 * delta_y else delta_y,
           width = if(sign(delta_x) < 0) -1 * delta_x else delta_x)  
    )
    
  }

  calc_y = function(ys, bt_df) {
    return(
      as.list(
        as.numeric(
          subset(bt_df, subset = bt_df$y == ys, select = "layer")[[1]]
          )
        )
      )
  }

  calc_data = function(t, bands, stars_obj) {
    cat(paste(Sys.time(), "; Time: ", if(is.na(t)) 1 else t, "...\n", sep = ""))
    
    bt_raster = if(!is.na(t)) as(if(!is.na(bands)) stars_obj[,,,bands,t, drop = TRUE] 
                                  else stars_obj[,,,t, drop = TRUE], "Raster") 
                else as(if(!is.na(bands)) stars_obj[,,,bands, drop = TRUE] 
                        else stars_obj[,,], "Raster")

    bt_df = as.data.frame(bt_raster, xy = TRUE)
    uy = as.list(unique(bt_df[, 2]))
    y_list = lapply(uy, calc_y, bt_df)
    
    return(y_list)
  }

  if(!is.na(tot_bands)) {
    length(json_out$raster_collection_tiles) = tot_bands
    for(bands in 1:tot_bands) {
      cat(paste(Sys.time(), "; Processing Band: ", bands, "...\n", sep = ""))
      json_out$raster_collection_tiles[[bands]]$extent = calc_extent(stars_obj, bands)

      times = as.numeric(dim(stars_obj)["time"])
      if(!is.na(times)) {
        t_start = seq(from = attr(stars_obj[,,,bands,], "dimensions")$time$offset, by = attr(stars_obj[,,,bands,], "dimensions")$time$delta, length.out = times)
        t_end = c(t_start[2:length(t_start)], t_start[length(t_start)] + attr(stars_obj[,,,bands,], "dimensions")$time$delta)
        json_out$raster_collection_tiles[[bands]]$start_times = as.list(as.character.POSIXt(t_start, format = "%Y-%m-%dT%T %Z"))
        json_out$raster_collection_tiles[[bands]]$end_times = as.list(as.character.POSIXt(t_end, format = "%Y-%m-%dT%T %Z"))
        data = lapply(as.list(1:times), calc_data, bands, stars_obj)
      } else {
        t_start = NA
        t_end = NA
        json_out$raster_collection_tiles[[bands]]$start_times = as.list(NA)
        json_out$raster_collection_tiles[[bands]]$end_times = as.list(NA)
        data = lapply(as.list(NA), calc_data, bands, stars_obj)
      }
      json_out$raster_collection_tiles[[bands]]$data = data
    }
  } else {
    cat(paste(Sys.time(), "; Processing Band: 1;\n", sep = ""))
    length(json_out$raster_collection_tiles) = 1
    json_out$raster_collection_tiles[[1]]$extent = calc_extent(stars_obj, NA)
    bands = NA
    times = as.numeric(dim(stars_obj)["time"])
    if(!is.na(times)) {
      t_start = seq(from = attr(stars_obj[,,,], "dimensions")$time$offset, by = attr(stars_obj[,,,], "dimensions")$time$delta, length.out = times)
      t_end = c(t_start[2:length(t_start)], t_start[length(t_start)] + attr(stars_obj[,,,], "dimensions")$time$delta)
      json_out$raster_collection_tiles[[1]]$start_times = as.list(as.character.POSIXt(t_start, format = "%Y-%m-%dT%T %Z"))
      json_out$raster_collection_tiles[[1]]$end_times = as.list(as.character.POSIXt(t_end, format = "%Y-%m-%dT%T %Z"))
      data = lapply(as.list(1:times), calc_data, bands, stars_obj)
    } else {
      t_start = NA
      t_end = NA
      json_out$raster_collection_tiles[[1]]$start_times = as.list(NA)
      json_out$raster_collection_tiles[[1]]$end_times = as.list(NA)
      data = lapply(as.list(NA), calc_data, bands, stars_obj)
    }
    json_out$raster_collection_tiles[[1]]$data = data
  }
  cat("\n")
  print(Sys.time())
  cat("Converted resulting stars object back to JSON!\n")

  return(json_out)
}

json2dim_mod = function(json_dim) {
  dim_num = NA
  if(json_dim == "band") dim_num = 3
  if(json_dim == "time") dim_num = 4
  
  return(dim_num)
}

#' @serializer unboxedJSON
#' @post /udf
run_UDF.json = function(req) {
  print(Sys.time())
  cat("Started executing at endpoint /udf\n")
  json_in = fromJSON(req$postBody, simplifyVector = FALSE)
  script_text = json2script(json_in)

  dim_mod = try(json2dim_mod(json_in$code$alt_dim), silent = T)
  if(class(dim_mod) == "try-error") dim_mod = 4 # Testing

  stars_in = json2stars(json_in)
  stars_out = run_script(stars_obj = stars_in, dim_mod = dim_mod, script_text = script_text)
  json_out = stars2json(stars_obj = stars_out, json_in = json_in)

  # Generate HTTP response for "backend" with body as the JSON in the file `json_out_file`
  print(Sys.time())
  cat("Generating response to HTTP request")
  
  return(json_out)
}

#' @serializer unboxedJSON
#' @post /udf/raw
run_UDF.json.raw = function(req) {
  json_in = jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
  script_text = json2script(json_in)

  stars_in = json2stars(json_in)
  stars_out = run_script_raw(stars_obj = stars_in, script_text = script_text)
  json_out = stars2json(stars_obj = stars_out, json_in = json_in)

  #Generate HTTP response for "backend"
  print(Sys.time())
  cat("Generating resposne to HTTP request")
  
  return(json_out)
}

# Florian's subsetting hack re-introduced
# this is just a temporal fix for an issue during subsetting stars objects with an variable that was defined not in the basenv environment
"[.stars" = function(x, i = TRUE, ..., drop = FALSE, crop = TRUE) {
  missing.i = missing(i)
  # special case:
  if (! missing.i && inherits(i, c("sf", "sfc", "bbox")))
    return(stars::st_crop(x, i, crop = crop))

  mc <- match.call(expand.dots = TRUE)
  # select list elements from x, based on i:
  d = attr(x, "dimensions")
  ed = stars:::expand_dimensions.dimensions(d)
  x = unclass(x)[i]
  # selects also on dimensions:
  if (length(mc) > 3) {
    mc[[1]] <- `[`
    if (! missing(i))
      mc[[3]] <- NULL # remove i
    mc[["drop"]] = FALSE
    for (i in names(x)) {
      mc[[2]] = as.name(i)
      x[[i]] = eval(mc, x, enclos = parent.frame())
    }
    mc0 = mc[1:3] # "[", x, first dim
    j = 3 # first dim
    for (i in names(d)) {
      mc0[[2]] = as.name(i)
      mc0[[3]] = mc[[j]]
      mc0[["values"]] = ed[[i]]
      d[[i]] = eval(mc0, d, enclos = parent.frame())
      j = j + 1
    }
  }
  if (drop)
    abind::adrop(stars::st_as_stars(x, dimensions = d))
  else
    stars::st_as_stars(x, dimensions = d)
}

# *************************************************************
# RESTful web service with data as a base64 encoded string ====
# *************************************************************

close_relevant_conn = function(con_description) {
  cno = as.numeric(
    rownames(
      as.data.frame(showConnections())[as.data.frame(showConnections())$description == con_description]
      )
    )
  for(c in cno) {
    con = try(getConnection(what = c), silent = TRUE)
    if(class(con) != "try-error") {
      close(con)
      cat("\nConnection(s) closed successfully!\n")
    } else {
      cat("\nNo connections with given description to close.\n")
      break
    }
  }
}

bin_unzip_string = function(string = "data/binary_udf/bin_data", file = TRUE) {
  cat("Decoding base64 encoded string...\n")
  
  if(file) base64enc::base64decode(file = string, output = file("temp.zip", "wb")) 
  else base64enc::base64decode(what = string, output = file("temp.zip", "wb"))
  
  while(!file.exists("temp.zip")) Sys.sleep(1)
  
  close_relevant_conn("temp.zip")
  cat("Finished decoding string; Starting to uncompress ZIP file...\n")
  utils::unzip(zipfile = "temp.zip", overwrite = T, exdir = "disk") # Works with Windows
  cat("Finished unzipping file; Removing ZIP file...\n")
  file.remove("temp.zip")
  cat("Finished deleting ZIP file\n")
}

bin_read_legend = function(legend) { #FLA: nothing returned?
  cat("Creating stars object...\n")
  timestamps = unique(legend$timestamp)
  bands = unique(legend$band)
  filewpaths = cbind(legend[,1], legend$filename)[,2]
  stars_obj = openEO.R.UDF::read_stars(filewpaths, along = list(band = bands, time = timestamps)) #FLA: read_stars needs to be present in this plumbed file!
  
  return(stars_obj)
}


#' @serializer unboxedJSON
#' @post /udf/binary
run_UDF.binary = function(req) {
  cat("Reading JSON...\n")
  post_body = jsonlite::fromJSON(req$postBody) # for use with plumber
  cat("Converted incoming JSON to R object\n")

  # bin_unzip_string(string = post_body$base64str, file = FALSE)
  bin_unzip_string(string = post_body$base64str, file = FALSE)

  cat("Reading legend...\n")
  legend = jsonlite::fromJSON(post_body$legend) #FLA: is this necessary? if the body as a whole is read from JSON the values
  # should be already a R list
  
  legend$timestamp = as.POSIXct(legend$timestamp)
  stars_in = bin_read_legend(legend)
  cat("Creating stars object from incoming data\n")
  unlink("disk", recursive = TRUE)
  cat("Deleted directory disk\n")

  script = post_body$code$code$source
  script = gsub("\"", "'", script)
  script = gsub("\r", "", script)
  
  cat("Applying UDF on incoming stars object...\n")
  stars_out = run_script_raw(stars_obj = stars_in, script_text = script)
  cat("Output stars object created\n")

  time_only = FALSE
  band_only = FALSE
  
  time_out = try(dim(stars_out)[["time"]], silent = TRUE)
  if(class(time_out) == "try-error") {
    time_out = 1
    time_only = TRUE
  }
  
  band_out = try(dim(stars_out)[["band"]], silent = TRUE)
  if(class(band_out) == "try-error") {
    band_out = 1
    band_only = TRUE
  }

  legend_out = matrix(ncol = ncol(legend), nrow = time_out * band_out)
  colnames(legend_out) = colnames(legend)
  legend_out = as.data.frame(legend_out)
  cat("Outgoing legend created\n")

  out_dir = "results"
  dir.create(out_dir)
  if(!time_only) time_vals = attr(stars_out, "dimensions")[["time"]]$values else time_vals = NA
  if(!band_only) band_vals = attr(stars_out, "dimensions")[["band"]]$values else band_vals = NA
  
  cat("Starting to write results...\n")
  for(time_num in 1:time_out)  {
    cat(paste("Time:", time_num, "\n", sep = " "))
    out_path = paste(out_dir, "/t_", time_num, sep = "")
    dir.create(out_path)
    for(band_num in 1:band_out) {
      cat(paste("Band:", band_num, "\n", sep = " "))
      filename = paste(out_path, "/b_", band_num, ".tif",  sep = "")
      
      if(!time_only && !band_only) {
        stars_subset = stars_out[,,,band_num, time_num, drop = T]
      } else if(time_only) {
        stars_subset = stars_out[,,,time_num, drop = T] 
      } else if(band_only) {
        stars_subset = stars_out[,,,band_num, drop = T]
      }
      
      stars::st_write(obj = stars_subset, dsn = filename)
      index = ((time_num - 1) * band_out) + band_num
      
      legend_out[index,] = c(index, 
                             filename, 
                             as.numeric(time_num), 
                             as.character.Date(time_vals[time_num]), 
                             as.numeric(band_num), 
                             band_vals[band_num])
    }
  }

  filepaths = list.files("results", full.names = T, recursive = T)
  zip::zip(zipfile = "results.zip", files = filepaths, recurse=TRUE)
  unlink("results", recursive = TRUE)
  
  out_bin_string = base64enc::base64encode(what = "results.zip")
  cat("Created outgoing base64 encoded string\n")
  
  file.remove("results.zip")
  
  response = list(legend = as.list(legend_out), base64str = out_bin_string)
  
  cat("Converted R object to JSON for response\n")

  return(response)
}


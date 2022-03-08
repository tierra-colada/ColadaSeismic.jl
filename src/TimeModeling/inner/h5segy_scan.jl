function segy_scan(files::Array{String,1}, keys::Array{String,1}, blocksize::Int; 
  chunksize::Int = CHUNKSIZE,
  pool::WorkerPool = WorkerPool(workers()),
  verbosity::Int = 1,
  filter::Bool = true)

  files_sort = files[sortperm(filesize.(files), rev = true)]
  run_scan(f) = scan_file(f, keys, blocksize, chunksize=chunksize, verbosity=verbosity)
  s = pmap(run_scan, pool, files_sort)

  return merge(s)
end
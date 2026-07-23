# Boots the Expression Console server: frontend + auth + analysis API on one port.
# Render/Cloud Run/App Runner inject $PORT; locally it falls back to 8000.
port <- as.integer(Sys.getenv("PORT", "8000"))
cat(sprintf(">>> Expression Console on 0.0.0.0:%d  (frontend + API)\n", port))

pr <- plumber::plumb("/app/plumber.R")
pr$run(host = "0.0.0.0", port = port, docs = FALSE)

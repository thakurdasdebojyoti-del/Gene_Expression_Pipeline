# Boots the plumber API. Render injects $PORT and expects the app to bind to it;
# locally it falls back to 8000. Binding to 0.0.0.0 is required on Render.
port <- as.integer(Sys.getenv("PORT", "8000"))
cat(sprintf(">>> Starting plumber API on 0.0.0.0:%d\n", port))

pr <- plumber::plumb("/app/plumber.R")
pr$run(host = "0.0.0.0", port = port, docs = FALSE)
